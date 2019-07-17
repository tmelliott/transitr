# 0.4 - fourth beta release

- now using NODES for the network, which can be stops or intersections
  (currently only stops are working, but hopefully flexible to extend later)
- major improvements to particle filter, now correctly estimating 
  travel times for all segments


# 0.3 - third beta release

- particle filter running better
- network running better
- trip ETAs generated per trip using KF (instead of PF)


# 0.2 - second beta release

- real-time vehicle model running with configurable parameters
- real-time network model "implemented" (but in a very basic stage)
- likewise, network model "implemented"
- verbosity turned off by default


# 0.1.0 - first beta release

- basic GTFS data import into a database
- connect to a live feed (with API key if required)
- model vehicles using a very simple model
- generate current-speed based ETAs
