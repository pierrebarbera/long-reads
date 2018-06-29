/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2018 Lucas Czech and HITS gGmbH

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
    Lucas Czech <lucas.czech@h-its.org>
    Exelixis Lab, Heidelberg Institute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "genesis/genesis.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <limits>
#include <string>
#include <unordered_map>

using namespace genesis;
using namespace genesis::placement;
using namespace genesis::tree;
using namespace genesis::utils;

int main( int argc, char** argv )
{
    // -------------------------------------------------------------------------
    //     Startup and command line args
    // -------------------------------------------------------------------------

    // Activate logging.
    utils::Logging::log_to_stdout();
    utils::Logging::details.time = true;
    LOG_BOLD << utils::Options::get().info();
    LOG_BOLD;

    LOG_INFO << "Started";

    // Check if the command line contains the right number of arguments.
    if (argc != 3) {
        throw std::runtime_error(
            "Need to provide two arguments: jplace file and output dir\n"
        );
    }

    // In out dirs.
    auto epafile = std::string( argv[1] );
    auto outdir = utils::dir_normalize_path( std::string( argv[2] ));
    utils::dir_create(outdir);

    // -------------------------------------------------------------------------
    //     Read sample file
    // -------------------------------------------------------------------------

    // Read data.
    LOG_INFO << "reading jplace sample";
    auto sample = JplaceReader().from_file( epafile );
    sort_placements_by_weight( sample );

    // -------------------------------------------------------------------------
    //     Preparation
    // -------------------------------------------------------------------------

    // Calc edpl.
    LOG_INFO << "edpl";
    auto const edpls = edpl( sample );
    auto const max_edpl = *std::max_element( edpls.begin(), edpls.end());
    LOG_INFO << "max edpl " << max_edpl;

    // Calc node dist mat.
    LOG_INFO << "dists";
    auto const node_distances = node_branch_length_distance_matrix( sample.tree() );
    auto const edge_distances = edge_path_length_matrix( sample.tree() );

    // Prepare result color vecs.
    auto const base_color = Color( 0.753, 0.753, 0.753 );
    auto colors_lwr  = std::vector<Color>( sample.tree().edge_count(), base_color );
    auto colors_edpl = std::vector<Color>( sample.tree().edge_count(), base_color );

    // Go crazy with color.
    // auto inferno = color_list_inferno();
    // inferno.erase( inferno.begin() + 4 * inferno.size() / 5, inferno.end() );
    auto inferno = std::vector<Color>{
        color_from_hex( "#00000C" ),
        color_from_hex( "#87216B" ),
        color_from_hex( "#C83E4F" ),
        color_from_hex( "#FFA02F" )
    };

    // Color for LWR
    auto map_lwr = ColorMap( inferno );
    // auto norm_lwr  = utils::ColorNormalizationLinear();
    // norm_lwr.scale( 0.0, 1.0 );
    auto norm_lwr  = utils::ColorNormalizationBoundary();
    norm_lwr.scale( 0.0, 1.0, 4 );
    map_lwr.reverse( true );

    // Color for EDPL
    auto map_edpl = ColorMap( inferno );
    // auto norm_edpl = utils::ColorNormalizationLinear();
    // norm_edpl.scale( 0.0, max_edpl );
    auto norm_edpl = utils::ColorNormalizationBoundary();
    norm_edpl.scale( 0.0, max_edpl, 4 );

    // count how many seqs did not have a placement on their branch.
    size_t no_booking = 0;

    // weighted distances to the branch where each seq is on the tree,
    // in branch length units (continuous) and as path lengths (discrete).
    auto weighted_bu_dists = std::vector<double>( sample.size(), 0.0 );
    auto weighted_pl_dists = std::vector<double>( sample.size(), 0.0 );

    // -------------------------------------------------------------------------
    //     Main Loop
    // -------------------------------------------------------------------------

    LOG_INFO << "queries";
    for( size_t qi = 0; qi < sample.size(); ++qi ) {
        auto const& pquery = sample.at( qi );

        // Find the branch that this pquery sits at.
        if( pquery.name_size() != 1 ) {
            LOG_ERR << "pquery.name_size() != 1";
            continue;
        }
        auto const& name = pquery.name_at(0).name;
        auto node_ptr = find_node( sample.tree(), name );
        if( ! node_ptr ) {
            LOG_ERR << "! node_ptr";
            continue;
        }
        auto const& node = *node_ptr;
        auto const br_idx = node.link().edge().index();

        // branch indices of the two branches adjacent to the query sequence branch,
        // as well as the branch itself (needed for safety checks only).
        auto const a_br_idx = node.link().edge().primary_link().next().edge().index();
        auto const b_br_idx = node.link().edge().primary_link().next().next().edge().index();
        auto const c_br_idx = node.link().edge().primary_link().next().next().next().edge().index();
        assert( node.link().edge().primary_link().next().next().next().next().edge().index() == a_br_idx );
        assert( node.link().edge().primary_link().edge().index() == c_br_idx );
        assert( node.link().edge().secondary_link().index() == node.link().index() );

        // node indices of the ends of the two branches next to the query branch.
        // this is what we measure distance to if the placement is not on one of those two branches.
        // in other words, we measure the distance of a placement to the beginning of one of the two
        // branches where the sequences is placed in the comprehensive tree.
        // this is a slight underestimation, as we ignore the positon of the sequence between the two
        // branches, but those are not reliable anyway, so this is the best we can do.
        auto const a_nd_idx = node.link().edge().primary_link().next().outer().node().index();
        auto const b_nd_idx = node.link().edge().primary_link().next().next().outer().node().index();
        auto const p_nd_idx = node.link().edge().primary_link().node().index();

        // LOG_DBG << name << ": br_idx(" << br_idx << "), a_br_idx(" << a_br_idx << "), b_br_idx(" << b_br_idx << "), a_nd_idx(" << a_nd_idx << "), b_nd_idx(" << b_nd_idx << "), p_nd_idx(" << p_nd_idx << ")";

        // Check all placements and measure their distance to the branch.
        double lwr = -1.0;
        double weighted_bu_dist = 0.0;
        double weighted_pl_dist = 0.0;
        for( auto const& place : pquery.placements() ) {

            // LWR and safty checks
            auto const lwrstr = to_string_precise( place.like_weight_ratio );
            if( place.edge().index() == a_br_idx ) {
                // LOG_DBG1 << "its a: " << lwrstr << " at " << name;
                if( lwr != -1.0 ) {
                    LOG_WARN << "double booking for " << name;
                }
                lwr = place.like_weight_ratio;
            }
            if( place.edge().index() == b_br_idx ) {
                // LOG_DBG1 << "its b: " << lwrstr << " at " << name;
                if( lwr != -1.0 ) {
                    LOG_WARN << "double booking for " << name;
                }
                lwr = place.like_weight_ratio;
            }
            if( place.edge().index() == c_br_idx ) {
                LOG_WARN << "impossible! its c: " << lwrstr << " at " << name;
            }
            if( place.edge().index() == br_idx ) {
                LOG_WARN << "impossible: " << lwrstr << " at " << name;
            }

            // LOG_DBG1 << lwrstr << " at edge index " << place.edge().index();

            // if the placement is on one of the good branches, we do not add any weight to the sum.
            if( place.edge().index() == a_br_idx || place.edge().index() == b_br_idx ) {
                continue;
            }

            auto const br_len = place.edge().data<PlacementEdgeData>().branch_length;

            // If not, add up min dist to one of the branches.
            double const ap_dist = place.proximal_length
                + node_distances( a_nd_idx, place.edge().primary_node().index() )
            ;
            double const ad_dist = br_len
                - place.proximal_length
                + node_distances( a_nd_idx, place.edge().secondary_node().index() )
            ;
            double const bp_dist = place.proximal_length
                + node_distances( b_nd_idx, place.edge().primary_node().index() )
            ;
            double const bd_dist = br_len
                - place.proximal_length
                + node_distances( b_nd_idx, place.edge().secondary_node().index() )
            ;

            // Because the leave one out test does not print proper distal lengths,
            // we also do the same thing just without using prox and dist lengths...
            // it is a bit hacky, but does it's job for now.
            // That is, we ignore the position of the placement on its branch,
            // and just measure the distance from the two end nodes of its branch
            // to the two nodes where the sequences is placed in the comprehensive tree.
            double const aph_dist = node_distances( a_nd_idx, place.edge().primary_node().index() );
            double const adh_dist = node_distances( a_nd_idx, place.edge().secondary_node().index() );
            double const bph_dist = node_distances( b_nd_idx, place.edge().primary_node().index() );
            double const bdh_dist = node_distances( b_nd_idx, place.edge().secondary_node().index() );

            // Finally: add up
            // auto const mind = std::min({ ap_dist, ad_dist, bp_dist, bd_dist });
            auto const mind = std::min({ ap_dist, ad_dist, bp_dist, bd_dist, aph_dist, adh_dist, bph_dist, bdh_dist });
            weighted_bu_dist += mind * place.like_weight_ratio;
            // LOG_DBG1 << mind;

            // Safty: measure dist to central node. should never be smaller.
            double const pp_dist = place.proximal_length
                + node_distances( p_nd_idx, place.edge().primary_node().index() )
            ;
            double const pd_dist = br_len
                - place.proximal_length
                + node_distances( p_nd_idx, place.edge().secondary_node().index() )
            ;
            if( std::min( pp_dist, pd_dist ) < mind ) {
                LOG_WARN << "std::min( pp_dist, pd_dist ) < mind";
            }

            // now do the same for discrete distances.
            double const apl_dist = edge_distances(
                place.edge().index(),
                a_br_idx
            );
            double const bpl_dist = edge_distances(
                place.edge().index(),
                b_br_idx
            );
            auto minpl = std::min( apl_dist, bpl_dist );
            weighted_pl_dist += minpl * place.like_weight_ratio;
        }

        // If we never set the lwr, the query has no placement mass on the good branches.
        if( lwr < 0.0 ) {
            ++no_booking;
            // LOG_WARN << "no booking for " << name;
        }
        // colors_lwr[ pquery.placement_at(0).edge().index() ] = Color( 1,0,0 );

        // Safety.
        if( colors_lwr[ br_idx ] != base_color ) {
            LOG_WARN << "colors_lwr[ br_idx ] != base_color";
        }
        if( colors_edpl[ br_idx ] != base_color ) {
            LOG_WARN << "colors_edpl[ br_idx ] != base_color";
        }

        // Store values at branch. cannot store the weighted dist as color yet,
        // because we need to wait for its max.
        weighted_bu_dists[qi] = weighted_bu_dist;
        weighted_pl_dists[qi] = weighted_pl_dist;
        colors_lwr[ br_idx ]  = map_lwr( norm_lwr, lwr );
        colors_edpl[ br_idx ] = map_edpl( norm_edpl, edpls[ qi ] );
    }

    // -------------------------------------------------------------------------
    //     Final calculations
    // -------------------------------------------------------------------------

    // report how many sequences have no placement on their branch
    LOG_INFO << no_booking << " out of " << sample.size() << " did not have a placement on their branch";

    // weighted bu dists: find max and define color map accordingly.
    auto const max_wd_bu = *std::max_element( weighted_bu_dists.begin(), weighted_bu_dists.end());
    auto map_wd_bu = ColorMap( inferno );
    // auto norm_wd_bu = utils::ColorNormalizationLinear();
    // norm_wd_bu.scale( 0.0, max_wd_bu );
    auto norm_wd_bu = utils::ColorNormalizationBoundary();
    norm_wd_bu.scale( 0.0, max_wd_bu, 4 );
    LOG_INFO << "max_wd_bu " << max_wd_bu;

    // weighted pl dists: find max and define color map accordingly.
    auto const max_wd_pl = *std::max_element( weighted_pl_dists.begin(), weighted_pl_dists.end());
    auto map_wd_pl = ColorMap( inferno );
    // auto norm_wd_pl = utils::ColorNormalizationLinear();
    // norm_wd_pl.scale( 0.0, max_wd_pl );
    auto norm_wd_pl = utils::ColorNormalizationBoundary();
    norm_wd_pl.scale( 0.0, max_wd_pl, 4 );
    LOG_INFO << "max_wd_pl " << max_wd_pl;

    auto colors_wd_bu  = std::vector<Color>( sample.tree().edge_count(), base_color );
    auto colors_wd_pl  = std::vector<Color>( sample.tree().edge_count(), base_color );

    // Fill weighted dist colors
    for( size_t qi = 0; qi < sample.size(); ++qi ) {
        auto const& pquery = sample.at( qi );
        auto const& name = pquery.name_at(0).name;
        auto node_ptr = find_node( sample.tree(), name );
        if( ! node_ptr ) {
            LOG_ERR << "! node_ptr";
            continue;
        }
        auto const& node = *node_ptr;
        auto const br_idx = node.link().edge().index();

        colors_wd_bu[ br_idx ] =  map_wd_bu( norm_wd_bu, weighted_bu_dists[qi] );
        colors_wd_pl[ br_idx ] =  map_wd_pl( norm_wd_pl, weighted_pl_dists[qi] );
    }

    // -------------------------------------------------------------------------
    //     Output files
    // -------------------------------------------------------------------------

    // Write color trees, circular.
    auto params = LayoutParameters();
    params.stroke.width = 8;
    write_color_tree_to_svg_file( sample.tree(), params, colors_lwr, map_lwr, norm_lwr, outdir + "lwr_c.svg" );
    write_color_tree_to_svg_file( sample.tree(), params, colors_edpl, map_edpl, norm_edpl, outdir + "edpl_c.svg" );
    write_color_tree_to_svg_file( sample.tree(), params, colors_wd_bu, map_wd_bu, norm_wd_bu, outdir + "wd_bu_c.svg" );
    write_color_tree_to_svg_file( sample.tree(), params, colors_wd_pl, map_wd_pl, norm_wd_pl, outdir + "wd_pl_c.svg" );

    // Rect
    params.shape = LayoutShape::kRectangular;
    params.type = LayoutType::kPhylogram;
    write_color_tree_to_svg_file( sample.tree(), params, colors_lwr, map_lwr, norm_lwr, outdir + "lwr_r.svg" );
    write_color_tree_to_svg_file( sample.tree(), params, colors_edpl, map_edpl, norm_edpl, outdir + "edpl_r.svg" );
    write_color_tree_to_svg_file( sample.tree(), params, colors_wd_bu, map_wd_bu, norm_wd_bu, outdir + "wd_bu_r.svg" );
    write_color_tree_to_svg_file( sample.tree(), params, colors_wd_pl, map_wd_pl, norm_wd_pl, outdir + "wd_pl_r.svg" );

    // Nexus
    write_color_tree_to_nexus_file( sample.tree(), colors_lwr, outdir + "lwr.nexus" );
    write_color_tree_to_nexus_file( sample.tree(), colors_edpl, outdir + "edpl.nexus" );
    write_color_tree_to_nexus_file( sample.tree(), colors_wd_bu, outdir + "wd_bu.nexus" );
    write_color_tree_to_nexus_file( sample.tree(), colors_wd_pl, outdir + "wd_pl.nexus" );

    LOG_INFO << "Finished";
    return 0;
}
