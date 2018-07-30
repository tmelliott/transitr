#include "realtime_feed.h"

static size_t WriteCallback(void *contents, size_t size, size_t nmemb, void *userp)
{
    ((std::string*)userp)->append((char*)contents, size * nmemb);
    return size * nmemb;
}

RealtimeFeed::RealtimeFeed (std::string& url, List& hdrs) : _url (url)
{
    _headers.reserve (hdrs.size ());
    for (int i=0; i<hdrs.size (); i++)
    {
        List hdr = hdrs[i];
        String hdrn = hdr["name"];
        String hdrv = hdr["value"];
        std::ostringstream hss;
        hss << (std::string)hdrn << ": " << (std::string)hdrv;
        _headers.push_back (hss.str ());
    }
}

int RealtimeFeed::update ()
{
    CURL *curl;
    CURLcode res;
    curl_global_init (CURL_GLOBAL_DEFAULT);
    curl = curl_easy_init ();

    _feed.Clear ();

    std::string readBuffer;
    if (curl)
    {
        // add the headers
        struct curl_slist *chunk = NULL;

        for (auto header: _headers)
        {
            chunk = curl_slist_append (chunk, header.c_str ());
        }

        res = curl_easy_setopt(curl, CURLOPT_HTTPHEADER, chunk);
        // curl_easy_setopt(curl, CURLOPT_VERBOSE, 1L);

        curl_easy_setopt (curl, CURLOPT_URL, _url.c_str ());
        curl_easy_setopt (curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt (curl, CURLOPT_WRITEDATA, &readBuffer);
        res = curl_easy_perform (curl);
        if (res != CURLE_OK)
        {
            Rcerr << "curl_easy_perform() failed\n";
            return 1;
            // curl_easy_strerror(res));
        }
        curl_easy_cleanup (curl);
        curl_slist_free_all (chunk);
    }
    curl_global_cleanup ();

    std::istringstream buf (readBuffer);
    if (!_feed.ParseFromIstream (&buf)) {
        Rcerr << "\n x Failed to parse GTFS realtime feed!\n";
        return 2;
    }

    return 0;
}

transit_realtime::FeedMessage* RealtimeFeed::feed ()
{
    return &_feed;
}

void load_vehicles (Gtfs::vehicle_map* vehicles,
                    transit_realtime::FeedMessage* feed,
                    Gtfs::Gtfs* gtfs)
{
    for (int i=0; i<feed->entity_size (); ++i)
    {
        auto ent = feed->entity (i);
        if (!ent.has_vehicle ()) continue;
        if (!ent.vehicle ().has_vehicle ()) continue;
        
        std::string id (ent.vehicle ().vehicle ().id ());
        auto vs = vehicles->find (id);
        if (vs == vehicles->end ())
        {
            auto r = vehicles->emplace (id, id);
            if (r.second)
            {
                r.first->second.update (ent.vehicle (), &(*gtfs));
            }
        }
        else
        {
            vs->second.update (ent.vehicle (), &(*gtfs));
        }
    }
}
