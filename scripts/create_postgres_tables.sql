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
)
