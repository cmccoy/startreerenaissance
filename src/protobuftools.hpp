#ifndef PROTOBUFTOOLS_H
#define PROTOBUFTOOLS_H

#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/io/gzip_stream.h>

#include <iosfwd>
#include <string>
#include <vector>

template<typename T>
std::vector<T> loadDelimitedFromStream(std::istream& in, const bool isGzipped=true)
{
    google::protobuf::io::IstreamInputStream raw_in(&in);
    google::protobuf::io::GzipInputStream zip_in(&raw_in);
    std::vector<T> result;
    google::protobuf::io::ZeroCopyInputStream* stream;
    if(isGzipped)
        stream = &zip_in;
    else
        stream = &raw_in;

    while(true) {
        google::protobuf::io::CodedInputStream coded_in(stream);
        uint32_t size = 0;
        bool success = coded_in.ReadVarint32(&size);
        if(!success) break;
        T item;
        std::string s;
        coded_in.ReadString(&s, size);
        success = item.ParseFromString(s);
        assert(success && "Failed to parse");
        result.push_back(std::move(item));
    }
    return result;
}

template<typename T>
void writeDelimitedToStream(std::ostream& out, const std::vector<T>& items, const bool gzip=true) {
    google::protobuf::io::OstreamOutputStream raw_out(&out);
    google::protobuf::io::GzipOutputStream zip_out(&raw_out);
    google::protobuf::io::CodedOutputStream coded_out(&zip_out);
    google::protobuf::io::ZeroCopyOutputStream* stream;
    if(gzip)
        stream = &zip_out;
    else
        stream = &raw_out;

    for(const T& item : items) {
        google::protobuf::io::CodedOutputStream coded_out(stream);
        coded_out.WriteVarint32(item.ByteSize());
        item.SerializeWithCachedSizes(&coded_out);
    }
}

#endif
