// Log.hh
// T. M. Kelley
// Jan 24, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

#include <queue>
#include <string>

namespace nut
{
    struct Std_Log
    {
        typedef std::string entry_t;

        void operator()(entry_t const & entry){ 
            std::cout << entry << std::endl;
        }

        bool isNull(){return false;}
    };

    struct Queue_Log
    {
        typedef std::string entry_t;
        typedef std::queue<entry_t> entries_t;

        entries_t entries;

        void operator()(entry_t const & entry){ entries.push(entry);}

        bool isNull(){return false;}
    };

    struct Null_Log
    {
        typedef std::string entry_t;
        typedef std::queue<entry_t> entries_t;

        void operator()(entry_t const &) {}

        bool isNull(){return true;}
    };

} // nut::

// version
// $Id$

// End of file
