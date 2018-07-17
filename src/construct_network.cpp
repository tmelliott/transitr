// [[Rcpp::plugins("cpp11")]]

#include <vector>
#include <Rcpp.h>

#include <vector>
#include <algorithm>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

/**
#include <sqlite3.h>

#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

// definitions
namespace bg = boost::geometry;
namespace bgi = bg::index;

// typedef bg::model::point<double, 2, bg::cs::cartesian> node_t;
typedef bg::model::point<double, 2, bg::cs::spherical_equatorial<bg::degree> > node_t;
typedef std::pair<node_t, unsigned int> value_t;

typedef bg::model::segment<node_t> edge_t;
typedef std::pair<edge_t, unsigned int> evalue_t;

#include "geo.h"
#include "node.h"

using namespace Rcpp;

int closest_node (Node& reference_node, std::vector<Node>* existing_nodes);
int add_node (std::vector<Node>* nodes, Node& node);
int add_edge (std::vector<Edge>* edges, int n1, int n2);
node_t project_point (node_t n, edge_t e);


// [[Rcpp::export]]
List construct_network (List nw, List shapes)
{
  // nw is used for writing the results to the database

  
  // NumericMatrix nodes;
  std::vector<Node> nodes;
  std::vector<Edge> edges;
  
  // Node Rtree
  bgi::rtree<value_t, bgi::quadratic<16> > node_tree;
  bgi::rtree<evalue_t, bgi::quadratic<16> > edge_tree;
  std::vector<unsigned> removed_edges;
  
  int prevNode;
  int curNode;
  int Nt (shapes.size ());
  int Ni (0);
  
  auto tstart = Clock::now ();
  for (auto& ti : shapes)
    {
      Ni++;
      Rcout << "\r ~ Processing shape " << Ni << " of " << Nt 
            << " - Network has " << node_tree.size () << " nodes and "
            << edge_tree.size () - removed_edges.size () << " edges";
      NumericMatrix shape = ti;
      prevNode = -1;
      for (unsigned int i=0; i<shape.nrow (); i++)
        {
	  // Rcout << "\n > Point " << (i+1) << ": ";
	  // Query Rtree for nodes within 10m of proposed nodes
	  std::vector<value_t> result_n;
	  node_t p (shape (i,0), shape (i,1));
	  node_tree.query (
			   bgi::nearest (p, 1) && 
			   bgi::satisfies ([&](value_t const& v) {
			       return bg::distance(v.first, p) * earthRadius < 5;
			     }), 
			   std::back_inserter (result_n));

	  if (result_n.size ()) {
	    // Found neaby node - merge new node with closest node
	    curNode = result_n[0].second;
	  } else {
	    // no matching nodes, look for nearby edges
	    // result is [matched edge, Projected point]
	    std::vector<std::tuple<evalue_t, node_t, double> > result_e;
	    for (auto it = edge_tree.qbegin (
					     bgi::satisfies ([&](evalue_t const& v) {
						 // find edge in removed_ids
						 if (std::find (removed_edges.begin (), removed_edges.end (), v.second)
						     != removed_edges.end ()) return false;
						 return bg::distance(v.first, p) * earthRadius < 10;
					       })); it != edge_tree.qend (); ++it)
	      {
		// get the projected point, and check that it lies between A and B

		// lon/lat
		std::pair<double,double> p0 (bg::get<0> (p), bg::get<1> (p));
		std::pair<double,double> p1 (bg::get<0,0> (it->first), bg::get<0,1> (it->first));
		std::pair<double,double> p2 (bg::get<1,0> (it->first), bg::get<1,1> (it->first));

		double length_e = bg::length (it->first) * earthRadius;
		double ctd = fabs(crossTrackDistance ( // lon, lat
						      std::get<1> (p0), std::get<0> (p0),
						      std::get<1> (p1), std::get<0> (p1),
						      std::get<1> (p2), std::get<0> (p2)));

		if (ctd < 0 || ctd > length_e) continue;

		// double check that projection distance is < 20
		double atd = alongTrackDistance (
						 std::get<1> (p0), std::get<0> (p0),
						 std::get<1> (p1), std::get<0> (p1),
						 std::get<1> (p2), std::get<0> (p2));
		double b (bearing (std::get<1> (p1), std::get<0> (p1), 
				   std::get<1> (p2), std::get<0> (p2)));
		Node P (destinationPoint (p1, b, atd));
		node_t Pn (std::get<0> (P), std::get<1> (P));

		double pdist (bg::distance (Pn, p) * earthRadius);
		if (pdist > 10) continue;

		// check that direction is the same
		// i.e., angle from prevNode to P - bearing(edge) (== b) < 45
		Node prevn (nodes[prevNode]);
		double b1 = bearing (std::get<1> (prevn), std::get<0> (prevn),
				     std::get<1> (P), std::get<0> (P));
		if (fabs (b - b1) > 45) continue;
                    
		result_e.emplace_back (*it, Pn, pdist);
	      }

	    if (result_e.size ()) {
	      // Rcout << "\n *** found " << result_e.size () << " candidate, using closest";
	      // merge with or split the closest edge
                    
	      // use the edge with the shortest projection distance
	      auto res = std::min_element (result_e.begin (), result_e.end (),
					   [] (std::tuple<evalue_t, node_t, double> s1, std::tuple<evalue_t, node_t, double> s2) {
					     return std::get<2> (s1) < std::get<2> (s2);
					   });
	      evalue_t e (std::get<0> (*res));
	      node_t n (std::get<1> (*res));
	      Node newnode (bg::get<0> (n), bg::get<1> (n));

	      // lon/lat
	      std::pair<double,double> p1 (bg::get<0,0> (e.first), bg::get<0,1> (e.first));
	      std::pair<double,double> p2 (bg::get<1,0> (e.first), bg::get<1,1> (e.first));
	      double d1 = distanceEarth(std::get<1> (p1), std::get<0> (p1),
					std::get<1> (newnode), std::get<0> (newnode));
	      double d2 = distanceEarth(std::get<1> (p2), std::get<0> (p2),
					std::get<1> (newnode), std::get<0> (newnode));

	      if (d1 < d2 && d1 < 5) {
		// merge P into p1
		curNode = std::get<0> (edges[e.second]);
	      } else if (d2 < d1 && d2 < 5) {
		// merge P into p2
		curNode = std::get<1> (edges[e.second]);
	      } else {
		curNode = add_node (&nodes, newnode);
		node_tree.insert (std::make_pair (n, static_cast<unsigned int> (curNode)));

		Edge e_old = edges[e.second];
		std::get<0> (edges[e.second]) = -1;
		std::get<1> (edges[e.second]) = -1;

		// split p1p2 into p1P and p2P
		edge_t e1 (node_t (std::get<1> (p1), std::get<0> (p1)), n);
		int edgeID = add_edge(&edges, std::get<0> (e_old), curNode);
		if (edgeID > -1)
		  edge_tree.insert (std::make_pair (e1, static_cast<unsigned int> (edgeID)));

		edge_t e2 (n, node_t (std::get<1> (p2), std::get<0> (p2)));
		edgeID = add_edge(&edges, curNode, std::get<1> (e_old));
		if (edgeID > -1)
		  edge_tree.insert (std::make_pair (e2, static_cast<unsigned int> (edgeID)));

		// delete p1p2
		removed_edges.push_back (e.second);
	      }

	    } else {
	      // no nearby edges, insert new node
	      // Rcout << "no matches - adding";
	      node_t n (shape (i,0), shape (i,1));
	      Node newnode (shape (i, 0), shape (i, 1));
	      curNode = add_node (&nodes, newnode);
	      node_tree.insert (std::make_pair (n, static_cast<unsigned int> (curNode)));
	    }
	  }
	  // Rcout << "\n";

	  if (prevNode > -1) {
	    // add edge if it doesn't exist
	    Node from = nodes[prevNode];
	    Node to = nodes[curNode];
	    node_t n1 (std::get<0> (from), std::get<1> (from));
	    node_t n2 (std::get<0> (to), std::get<1> (to));
	    edge_t e (n1, n2);
	    int edgeID = add_edge(&edges, prevNode, curNode);
	    if (edgeID > -1)
	      edge_tree.insert (std::make_pair (e, static_cast<unsigned int> (edgeID)));
	  }
	  prevNode = curNode;
        }
    }

  auto tend = Clock::now ();
  Rcout << "\r Processing " << Nt << " shapes took "
        << std::chrono::duration<double, std::milli> (tend - tstart).count () / 1000
        << " seconds             \n"
        << " * " << node_tree.size () << " nodes\n"
        << " * " << edge_tree.size () - removed_edges.size () << " edges\n";

  NumericMatrix Mnodes (nodes.size (), 2);
  for (int i=0; i<nodes.size (); i++)
    {
      NumericMatrix::Row ri = Mnodes ( i , _ );
      ri[0] = std::get<0> (nodes[i]);
      ri[1] = std::get<1> (nodes[i]);
    }

  IntegerMatrix Medges (edge_tree.size () - removed_edges.size (), 2);
  int i=0;
  for (auto it = edge_tree.qbegin ( bgi::satisfies ([&](evalue_t const& v) { return true; })); 
       it != edge_tree.qend (); ++it)
    {
      if (std::find (removed_edges.begin (), removed_edges.end (), it->second)
	  != removed_edges.end ()) continue;
      IntegerMatrix::Row ri = Medges ( i++ , _ );
      Edge e = edges[it->second];
      ri[0] = std::get<0> (e) + 1;
      ri[1] = std::get<1> (e) + 1;
    }

  List res = List::create();
  res.push_back (Mnodes, "nodes");
  res.push_back (Medges, "edges");
  res.attr("class") = "network";

  return res;
}

int add_node (std::vector<Node>* nodes, Node& node)
{
  nodes->push_back (node);
  return nodes->size () - 1;
}

int add_edge (std::vector<Edge>* edges, int n1, int n2)
{
  bool pathexists = false;
  int steps = 5;
  std::vector<int> to;
  to.push_back (n2);

  int M = edges->size ();

  // Rcout << "Looking for edge [" << n1 << ", " << n2 << "] ";
  while (!pathexists && steps > 0)
    {
      // collecting all possibilities
      // Rcout << "\n - Looking for edges terminating at nodes: ";
      // for (auto t : to) Rcout << t << ",";
      std::vector<int> candidates;
      for (int i=0; i<M; i++)
        {
	  Edge ei = (*edges)[i];
	  for (int j=0; j<to.size (); j++)
            {
	      if (std::get<1> (ei) == to[j])
                {
		  // Rcout << "\n   - Edge " << i << " [" 
		  // << std::get<0> (ei) << " -> " << std::get<1> (ei) << "]";
		  if (std::get<0> (ei) == n1) {
		    // Rcout << " MATCH";
		    pathexists = true;
		    break;
		  }
		  candidates.push_back (i);
                }
            }
        }
      if (candidates.size () == 0) break;
      to = candidates;
      steps--;
    }

  if (!pathexists) {
    // Rcout << "\n ** adding edge [" << n1 << ", " << n2 << "]\n";
    edges->emplace_back (n1, n2);
    return edges->size () - 1;
  }
  return -1;
}

node_t project_point (node_t n, edge_t e)
{
  return n;
}

**/