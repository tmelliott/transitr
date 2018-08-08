DROP TABLE IF EXISTS vehicles;
CREATE TABLE vehicles (
    vehicle_id VARCHAR(30) PRIMARY KEY NOT NULL,
    trip_id VARCHAR(80),
    timestamp TIMESTAMP,
    position_latitude REAL,
    position_longitude REAL,
    distance REAL,
    speed REAL,
    progress SMALLINT
);
DROP TABLE IF EXISTS vehicles_history;
CREATE TABLE vehicles_history (
    vehicle_id VARCHAR(30) NOT NULL,
    trip_id VARCHAR(80),
    timestamp TIMESTAMP,
    position_latitude REAL,
    position_longitude REAL,
    distance REAL,
    speed REAL,
    progress SMALLINT,
    PRIMARY KEY (vehicle_id, timestamp)
);
