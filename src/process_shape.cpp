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

    NumericMatrix tmp (L, 2);
    int k (0);
    tmp (k, _ ) = x (0, _ );

    double d;
    for (int i=1; i<L; i++)
    {
        d = distanceEarth(tmp (k, 1), tmp (k, 0), x (i, 1), x (i, 0));
        if (d  < 5) continue;
        k++;
        tmp (k, _ ) = x (i, _ );
    }

    NumericMatrix out (k + 1, 2);
    out = tmp ( Range (0, k), _ );

    out.attr ("class") = "network.shape";
    out.attr ("id") = id;
    return out;
}

