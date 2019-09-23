#ifndef FILE_READER_H
#define FILE_READER_H

#include "lattice_accessor.h"

#include <cstdio>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <unistd.h>
#include <memory>
#include <cassert>
#include <cstdlib>
#include <regex>
#include <map>

/*  
 *  These follow a bridge pattern with prototype IReader, concrete classes filetype_reader and bridge Reader.
 *   
 *  Writing a new parser?
 *      - Add your Readable_filetype to the Readable_filetype enumerator and extension_map below (Readable_file::, private)
 *      - Write your parser and (publicly) inherit and implement IReader
 *      - Add your class to Reader::read_objects_in()
 *      - Good to go!
 */ 

enum class Readable_filetype {
            NONE,
            VTK_STRUCTURED_GRID,
            PRO
        };

//Provides necessary checks before reading file by reader.
class Readable_file {
    public:
        const std::string m_filename;
        
        enum error {
            ERROR_EXTENSION,
            ERROR_FILE_NOT_FOUND,
        };

        Readable_filetype get_filetype();
        Readable_file(const std::string filename_, Readable_filetype filetype_);

    private:

        static std::map<Readable_filetype, std::string> extension_map;

        Readable_filetype m_filetype;
        std::string m_extension;

        void check_filetype();
        void read_extension();
};

//Interface for different filetypes
class IReader {
    public:
        IReader(Readable_file file);
        virtual ~IReader();

        virtual std::vector<std::vector<double>> get_file_as_vectors() = 0;
        std::vector<std::string> m_headers;

    protected:
        std::ifstream m_file;
        Lattice_accessor file_lattice;
        enum error {
            ERROR_FILE_FORMAT
        };

        virtual void set_lattice_geometry(const std::vector<std::string>&) = 0;

        virtual std::vector<std::string> tokenize(std::string line, char delimiter);

    private:
        //Disable default constructor
        IReader();
};

class Pro_reader : public IReader {
    private:
       std::vector<std::vector<double>> m_data;

        void check_delimiter(const std::string& line);
        void read_dimensions(const std::vector<std::string>& header_tokens);
        void check_component_name_format(const std::string& header_token);
        std::vector<std::string> parse_data(const size_t number_of_components, const size_t first_component_column);
        void set_lattice_geometry(const std::vector<std::string>& last_line);
        void adjust_indexing();

    public:
        explicit Pro_reader(Readable_file file);
        
        std::vector<std::vector<double>> get_file_as_vectors();
};


class Vtk_structured_grid_reader : public IReader {
    private:
        enum STATUS {
            END,
            NEW_BLOCK_FOUND,
            ERROR
        };

        void set_lattice_geometry(const std::vector<std::string>& tokens);
        STATUS parse_next_data_block(std::vector<double>& data);
        std::vector<double> with_bounds(std::vector<double>& input);

    public:

        Vtk_structured_grid_reader(Readable_file file);
        std::vector< std::vector<double> > get_file_as_vectors();
};

class Reader {
    public:
        Reader();
        ~Reader();

        //Callable for multiple files. returns number of objects read.
        size_t read_objects_in(Readable_file file);
        void push_data_to_objects(std::vector< std::vector<double> >& output);
        std::vector<std::string> get_headers() {return m_input_reader->m_headers;};


    private:
        std::vector< std::vector<double> > m_read_objects;
        std::ifstream m_file;
        std::unique_ptr<IReader> m_input_reader;
};


#endif