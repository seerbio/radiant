//
// Created by anichols on 8/23/23.
//

#include "Isoform.h"

Isoform::Isoform(const std::string &id): id(id) {}

Isoform::Isoform(
        const std::string &id,
        const std::string &name,
        const std::string &gene,
        const std::string &description,
        bool swissprot
        )
: id(id)
, name(name)
, gene(gene)
, description(description)
{}

