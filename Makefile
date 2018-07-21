proto: src/gtfs-realtime.pb.h

src/gtfs-realtime.pb.h: src/vendor/protobuf/gtfs-realtime.proto
	@protoc -I=src --cpp_out=./src/ $<
	@mv src/vendor/protobuf/*.cc src/
	