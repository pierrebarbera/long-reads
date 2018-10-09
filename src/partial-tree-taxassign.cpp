/*
    Copyright (C) 2018 Pierre Barbera

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Pierre Barbera <pierre.barbera@h-its.org>
    Exelixis Lab, Heidelberg Institute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "genesis/genesis.hpp"

#include <string>

using namespace genesis;
using namespace genesis::sequence;
using namespace genesis::tree;
using namespace genesis::utils;
using namespace genesis::taxonomy;

constexpr char UNDETERMINED[]   = "N/A";
constexpr char QUERY[]          = "query";

void print_labelled(Tree const& tree,
                    std::vector<Taxopath> const& node_labels )
{
    DefaultTreeNewickWriter writer;
    writer.node_to_element_plugins.push_back(
        [&]( TreeNode const& node, NewickBrokerElement& element ){
            element.comments.emplace_back(
                TaxopathGenerator().to_string( node_labels[node.index()] )
            );
        }
    );
    writer.to_stream( tree, std::cout );
}

Taxopath intersect( Taxopath const& lhs, Taxopath const& rhs )
{
    Taxopath result;

    // short-circuit if one is coming from a label: in that case ignore the query label
    if (lhs.size() and lhs[0] == QUERY) {
        return rhs;
    }
    if (rhs.size() and rhs[0] == QUERY) {
        return lhs;
    }

    // normal assignment if not
    for (size_t i = 0; ( i < std::min( lhs.size(), rhs.size() ) ) and ( lhs[i] == rhs[i] ); ++i) {
        result.push_back( lhs[i] );
    }

    if ( result.empty() ) {
        result.push_back( UNDETERMINED );
    }

    return result;
}

std::vector<Taxopath> label_nodes(  Tree const& tree,
                                    std::string const& taxon_file)
{
    TaxopathParser tpp;
    CsvReader csv_reader;
    csv_reader.separator_chars( "\t" );
    std::vector<Taxopath> node_labels(tree.node_count(), Taxopath({QUERY}));

    utils::InputStream it( utils::make_unique< utils::FileInputSource >( taxon_file ));
    while (it) {
        auto fields = csv_reader.parse_line( it );

        if ( fields.size() != 2 ) {
            throw std::runtime_error{"A line in the taxon file didn't have two tab separated columns."};
        }

        auto name = fields[0];
        std::string tax_string = fields[1];

        auto node_ptr = find_node( tree, name );

        if ( node_ptr == nullptr ) {
            throw std::runtime_error{"Could not find node with name: " + name};
        }

        node_labels[ node_ptr->index() ] = tpp.from_string( tax_string );
    }

    // check if any leafs weren't assigned a Taxopath
    // for ( auto const& node_it : tree.nodes() ) {
    //     if ( node_it->is_leaf() and node_labels[ node_it->index() ].empty() ) {
    //         auto name = node_it->data< DefaultNodeData >().name;
    //         throw std::runtime_error{"The leaf in the tree labelled '" + name
    //             + "' wasn't assigned a taxonomic path. Did you forget to include it in the taxon file?"};
    //     }
    // }
    // go through the tree in postorder fashion and label inner nodes according to the most common taxonomic rank of the children
    for ( auto it : postorder(tree) ) {
        if ( it.node().is_inner() ) {
            auto const child_1_idx = it.node().link().next().outer().node().index();
            auto const child_2_idx = it.node().link().next().next().outer().node().index();

            assert( not node_labels[ child_1_idx ].empty() );

            node_labels[ it.node().index() ] = intersect( node_labels[ child_1_idx ], node_labels[ child_2_idx ] );
        }
    }

    return node_labels;
}

void print_query_taxassign( std::ostream& stream,
                            Tree const& tree,
                            std::vector<Taxopath>& node_labels)
{
    // find a tip node that is labeled as the root point of the traversal
    TreeNode const * tmp_root = nullptr;
    // also get all the query tip node indices for later already
    std::vector<size_t> query_tip_indices;
    for ( auto const i : leaf_node_indices( tree ) ) {
        if ( node_labels[i][0] != QUERY ) {
            if (not tmp_root) {
                tmp_root = &tree.node_at( i );
            }
        } else {
            query_tip_indices.push_back( i );
        }
    }

    if ( tmp_root == nullptr ) {
        throw std::runtime_error{"None of the tips had labels!"};
    }

    // label all inner nodes via pre-order traversal
    TreeNode const * previous = tmp_root; // should work in our case?
    for ( auto it : tree::preorder(*tmp_root) ) {
        size_t i = it.node().index();
        if ( node_labels[i][0] == QUERY ) {
            // need to label according to parent!
            node_labels[i] = node_labels[ previous->index() ];
        }

        previous = &it.node();
    }

    // print all the query labels
    for ( auto const i : query_tip_indices ) {
        // output sativa-style taxassign
        stream << tree.node_at(i).data<DefaultNodeData>().name;
        stream << "\t" << TaxopathGenerator().to_string( node_labels[i] );
        // stream << "\t" << join( confidences, ";" );
        stream << "\n";
    }
}

/**
 * Takes a tree and a taxonomy file, which does only label a subset of the Trees taxa taxonomically.
 * Applies taxonomic labelling based on this partially labelled tree to the unlabelled queries, and prints
 * this information in tab-separated form to stdout
 */
int main( int argc, char** argv )
{
    // Check if the command line contains the right number of arguments.
    if (argc != 3) {
        throw std::runtime_error(
            std::string( "Usage: " ) + argv[0] + " <tree_file> <taxonomy_file>"
        );
    }

    std::string tree_file(argv[1]);
    std::string taxon_file(argv[2]);

    auto tree = DefaultTreeNewickReader().from_file( tree_file );

    auto node_labels = label_nodes(tree, taxon_file);

    print_query_taxassign( std::cout, tree, node_labels );

    // print_labelled(tree, node_labels);

    return 0;
}
