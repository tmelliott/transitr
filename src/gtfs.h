#ifndef GTFS_H
#define GTFS_H

#include <Rcpp.h>
#include <string>

class Gtfs {
private:
    std::string _dbname;

public:
    Gtfs (std::string& name) : _dbname (name) {
        Rcpp::Rcout << "Connected to GTFS database `"
            << _dbname << "`" << std::endl;
    };

    std::string& dbname () {
        return _dbname;
    }
};

#endif