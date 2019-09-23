#include "expander.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>

using namespace boost;
using namespace boost::program_options;

constexpr uint8_t BOUNDARIES = 2;

#include <iostream>
#include <fstream>
#include <exception>
#include <algorithm>
#include <iterator>
#include <vector>
#include <numeric>
#include <iomanip>
#include <random>

using namespace std;

int main(int argc , char **argv)
{
    options_description desc("\nExpand one-dimensional .pro files to multiple dimensions.\nFlips x-gradient to z.\nAllowed arguments");

    desc.add_options()
        ("help,h", "Print this help text.")
        ("dimensions,d", value< int >()->default_value(2), "[int] Number of dimensions to expand to.")
        ("y-dimension,y", value< int >(), "[int] Y size to expand to (without bounds).")
        ("z-dimension,z", value< int >(), "[int] Z size to expand to (without bounds).")
        ("preserve,p", value< bool >()->default_value(false), "[bool] Preserve bounds in the dimension being read.")
        ("input-file,i", value< string >(), "Specifies input file, must be pro format.")
        ("out-type,o", value< string >()->default_value("vtk"), "Specifies output file type (vtk_structured_grid, vtk_structured_points, or pro).")
        ("theta,t", value<vector<size_t>>()->multitoken(), "[int] [int] ... Sum and print theta of multiple components: index starts at 0, separated by spaces. Print multiple theta's by using this flag multiple times, e.g. -t 0 1 -t 2.")
        ("noise,n", value< double >(), "[double] Adds noise of given stddev to your perfectly smooth equilibrium profiles.");

    // Map positional parameters to their tag valued types 
    positional_options_description p;
    p.add("input-file", -1);

    Lattice_accessor lattice;

    // Parse the command line catching and displaying any 
    // parser errors
    variables_map vm;
    try
    {
        store( command_line_parser( argc, argv).options(desc).positional(p).run(), vm );
        notify(vm);

    } catch (std::exception &e)
    {
        cerr << endl << e.what() << endl;
        cerr << desc << endl;
    }

    if (vm.count("help")) {
        cerr << desc << endl;
        exit(0);
    }

    if (!vm.count("input-file")) {
        cerr << "No input file specified." << endl;
        exit(0);
    }

    if (!vm.count("y-dimension") or !vm.count("z-dimension")) {
        cerr << "Please specify both a y dimension size (-y arg) and z dimension size (-z arg). Bounds will be added automatically." << endl;
        exit(0);
    } else {
        if (vm["preserve"].as< bool >() == true) {
            lattice.MY = vm["y-dimension"].as< int >()-2;
            lattice.MZ = vm["z-dimension"].as< int >()-2;
        } else {
            lattice.MY = vm["y-dimension"].as< int >();
            lattice.MZ = vm["z-dimension"].as< int >();
        }
        lattice.dimensionality = static_cast<Dimensionality>(3);
    }

    auto map_it = Profile_writer::output_options.find(vm["out-type"].as< string >());

    if (map_it == Profile_writer::output_options.end())
        cerr << "Output type not recognized, please refer to help file" << endl;
        
    filesystem::path filename = vm["input-file"].as< string >();

    //else
    Readable_file in_file(filename.string(),Readable_filetype::PRO);

    Reader* in_reader = new Reader;
    size_t num_read_objects = in_reader->read_objects_in(in_file);

    // Should really have used a deque here, but namics only likes vectors so this
    // is compatible with the file reader from namics
    vector<vector<double>> input_densities(num_read_objects);
    in_reader->push_data_to_objects(input_densities);

    vector<string> headers = in_reader->get_headers();

    delete in_reader;

    /***** MAKE MIRROR IMAGE OF 1D SYSTEM IN X DIRECTION *****/
    if (vm["dimensions"].as< int >() == 3) {

        for (vector<double>& density_vector : input_densities) {
          if (vm["preserve"].as< bool >() == false) { 
            //remove boundary at position 0
            density_vector.erase(density_vector.begin());
          }

          vector<double> original = density_vector;

          //mirror densities
          std::reverse(density_vector.begin(), density_vector.end());
          density_vector.insert(density_vector.end(), original.begin(), original.end());
        }    
    }

    parsed_options parsed_options = command_line_parser(argc, argv)
        .options(desc)
        .run();

    std::vector<std::vector<std::string>> lists;
    for (const option& o : parsed_options.options) {
        if (o.string_key == "theta")
            lists.push_back( o.value );
    }

    variables_map vml;
    store(parsed_options, vml);

    for (auto list : lists) {
        double sum_theta_output{0.0};

        for (auto component_string : list) {
            size_t component = atoi(component_string.c_str());
            if (component < input_densities.size()) {
                sum_theta_output = std::accumulate(input_densities[component].begin(), input_densities[component].end(), sum_theta_output);
            } else {
                cerr << "Component " << component << " is out of range! Exiting." << endl;
                exit(0);
            }
        }

        if ( vm["preserve"].as< bool >() == true ) {
            sum_theta_output *= (lattice.MZ+BOUNDARIES)*(lattice.MY+BOUNDARIES);
        } else {
            sum_theta_output *= (lattice.MZ)*(lattice.MY);
        }
        cout << setprecision(14) << "(Sum) theta for component(s) ";
        for (auto component : list) {
            cout << component << " : ";
        }

            cout << sum_theta_output << endl;
    }

    // Used later on for writing, MY and MZ are set above
    lattice.MX = input_densities.back().size()-BOUNDARIES;
    lattice.set_jumps();

    /***** MULTIPLY SYSTEM IN Y AND Z DIRECTIONS *****/
    vector<vector<double>> output_densities(input_densities.size());
    int stride = (lattice.MZ+BOUNDARIES)*(lattice.MY+BOUNDARIES);

    for (size_t j = 0 ; j < input_densities.size() ; ++j) {
        output_densities[j].resize( lattice.system_size );
        for (size_t z = 0 ; z < input_densities[j].size() ; ++z) {
            for (size_t i = z*stride ; i < z*stride+stride ; ++i) {
                output_densities[j][i] = input_densities[j][z];
            }
        }
    } 

    if (vm.count("noise")) {
        /***** NOISE THE CRAP OUT OF THE PROFILE *****/
        std::minstd_rand prng { std::random_device{}() };
        double mean = 0;
        double stddev =  vm["noise"].as< double >();
        std::normal_distribution<double> dist{mean, stddev};

        for (auto& profile : output_densities) {
            for (auto& value : profile) {
                int output_index = (rand() % static_cast<int>(profile.size() + 1));
                double shift_quantity = value * dist(prng);
                double output_quantity = profile[output_index] + shift_quantity;
                double input_quantity = value - shift_quantity;

                while (not (output_quantity < 1.0 or output_quantity > 0.0) and (input_quantity < 1.0 or input_quantity > 0.0)) {
                    output_index = (rand() % static_cast<int>(profile.size() + 1));
                    shift_quantity = value * dist(prng);
                    output_quantity = profile[output_index] + shift_quantity;
                    input_quantity = value - shift_quantity;
                }

                profile[output_index] += shift_quantity;
                value -= shift_quantity;
            }
        }
    }


/* 
     NAMICS IS NOT ROW MAJOR BUT DEPTH MAJOR, SO THE GRADIENT WILL NOW BE IN Z INSTEAD OF X 

    size_t stride = (lattice.MX+BOUNDARIES)*(lattice.MY+BOUNDARIES);
    vector<vector<double>> output_densities(input_densities.size());

    for (size_t j = 0 ; j < input_densities.size() ; ++j) {
        output_densities[j].resize( lattice.system_size );
        for (size_t n = 0 ; n < stride; ++n)
            for (size_t i = 0 ; i < input_densities[j].size() ; ++i)
                output_densities[j][n*input_densities[j].size()+i] = input_densities[j][i];
    } */


    /***** PREPARE OUTPUT *****/
    auto out_filetype = Profile_writer::output_options[  vm["out-type"].as< string >()];

    Writable_file out_file(filename.stem().string() + "_expanded", out_filetype);
    auto profile_writer = Profile_writer::Factory::Create(out_filetype, &lattice, out_file);

    if (vm["preserve"].as< bool >() == true) {
        profile_writer->configuration.boundary_mode = IProfile_writer::Boundary_mode::WITH_BOUNDS;
    }

    for (size_t i = 0 ; i < output_densities.size() ; ++i)
        register_output_profile(headers[i], output_densities[i].data());

    profile_writer->bind_data(output_profiles);

    /***** WRITE HEADERS & DATA *****/
    profile_writer->prepare_for_data();

    profile_writer->write();

}