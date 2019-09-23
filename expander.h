#ifndef EXPANDER_H
#define EXPANDER_H

#include "file_reader.h"
#include "file_writer.h"
#include "lattice_accessor.h"

#include <memory>
#include <map>
#include <string>

std::map<std::string, std::shared_ptr<IOutput_ptr> > output_profiles;

template <typename Datatype>
void register_output_profile(std::string description, Datatype* variable) {
    std::shared_ptr<IOutput_ptr> profile = std::make_shared<Output_ptr<Datatype>>(variable);
    output_profiles[description] = profile;
}

#endif