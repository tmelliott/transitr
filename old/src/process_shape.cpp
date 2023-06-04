#include <Rcpp.h>

#include "geo.h"

using namespace Rcpp;

NumericMatrix process_shape (NumericMatrix x, String id);

// [[Rcpp::export]]
List shapes_df_to_list (DataFrame x)
{
    IntegerVector id = x[0];
    NumericVector lat = x[1];
    NumericVector lng = x[2];
    NumericVector seq = x[3];

    CharacterVector ids = id.attr ("levels");
    int N = ids.size ();
    int M = id.size ();

    List out (N);
    for (int j=0; j<N; j++)
    {
        std::vector<double> latj, lngj, seqj;
        latj.reserve (M);
        lngj.reserve (M);
        seqj.reserve (M);
        for (int i=0; i<M; i++)
        {
            if (id[i] != j+1) continue;
            latj.push_back (lat[i]);
            lngj.push_back (lng[i]);
            seqj.push_back (seq[i]);
        }
        NumericVector ylat = wrap (latj);
        NumericVector ylng = wrap (lngj);
        NumericVector yseq = wrap (seqj);
        latj.clear ();
        lngj.clear ();
        seqj.clear ();

        int L = ylat.size ();
        NumericMatrix y (L, 2);
        for (int i=0; i<L; i++)
        {
            // place them in the order of seq(uence)
            y (yseq[i]-1, 0) = ylng[i];
            y (yseq[i]-1, 1) = ylat[i];
        }

        out[j] = process_shape (y, ids[j]);
    }

    out.attr ("class") = "network.shape.list";

    return (out);
}

NumericMatrix process_shape (NumericMatrix x, String id)
{
    int L (x.nrow ());

    // NumericMatrix out (L, 3);
    NumericVector d (L);

    // out (0, Range (0, 1) ) = x (0, _ );
    // out (0, 2) = d;

    d[0] = 0.0;
    for (int i=1; i<L; i++)
    {
        d[i] = d[i-1] + 
            distanceEarth(x (i-1, 1), x (i-1, 0), x (i, 1), x (i, 0));
        // out (i, Range (0, 1) ) = x (i, _ );
        // out (i, 2) = d;
    }

    // bind distance to coordinates
    NumericMatrix out (L, 3);
    out (_, 0) = x (_,0);
    out (_, 1) = x (_,1);
    out (_, 2) = d;

    out.attr ("class") = "network.shape";
    out.attr ("id") = id;
    return out;
}

// [[Rcpp::export]]
NumericVector calculate_shape_distance (NumericMatrix x)
{
    int L (x.nrow ());
    // NumericVector lat = x( _ , 0);
    // NumericVector lng = x( _ , 1);
    // NumericVector seq = x ( _ , 2);
    // NumericMatrix y (L, 2);
    // for (int i=0; i<L; i++)
    // {
    //     // place them in the order of seq(uence)
    //     y (seq[i]-1, 0) = lat[i];
    //     y (seq[i]-1, 1) = lng[i];
    // }

    NumericVector d (L);
    d[0] = 0.0;
    for (int i=1; i<L; i++)
    {
        d[i] = d[i-1] + 
            distanceEarth(x (i-1, 1), x (i-1, 0), x (i, 1), x (i, 0));
    }

    return d;
}
