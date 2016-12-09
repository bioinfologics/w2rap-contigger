//
// Created by Bernardo Clavijo (TGAC) on 07/07/2016.
//

#include <paths/long/large/Lines.h>
#include "GFADump.h"

void GFADump (std::string filename, const HyperBasevector &hb, const vec<int> &inv, const
ReadPathVec &paths, const int MAX_CELL_PATHS, const int MAX_DEPTH, bool find_lines, const std::vector<int> & marked_edges){

    std::vector<std::string> colour_names={
            "aliceblue",
            "antiquewhite",
            "aqua",
            "aquamarine",
            "azure",
            "beige",
            "bisque",
            "blanchedalmond",
            "blue",
            "blueviolet",
            "brown",
            "burlywood",
            "cadetblue",
            "chartreuse",
            "chocolate",
            "coral",
            "cornflowerblue",
            "cornsilk",
            "crimson",
            "cyan",
            "darkblue",
            "darkcyan",
            "darkgoldenrod",
            "darkgreen",
            "darkgrey",
            "darkkhaki",
            "darkmagenta",
            "darkolivegreen",
            "darkorange",
            "darkorchid",
            "darkred",
            "darksalmon",
            "darkseagreen",
            "darkslateblue",
            "darkslategrey",
            "darkturquoise",
            "darkviolet",
            "deeppink",
            "deepskyblue",
            "dimgrey",
            "dodgerblue",
            "firebrick",
            "floralwhite",
            "forestgreen",
            "fuchsia",
            "gainsboro",
            "ghostwhite",
            "gold",
            "goldenrod",
            "grey",
            "green",
            "greenyellow",
            "honeydew",
            "hotpink",
            "indianred",
            "indigo",
            "ivory",
            "khaki",
            "lavender",
            "lavenderblush",
            "lawngreen",
            "lemonchiffon",
            "lightblue",
            "lightcoral",
            "lightcyan",
            "lightgoldenrodyellow",
            "lightgreen",
            "lightgrey",
            "lightpink",
            "lightsalmon",
            "lightseagreen",
            "lightskyblue",
            "lightslategrey",
            "lightsteelblue",
            "lightyellow",
            "lime",
            "limegreen",
            "linen",
            "magenta",
            "maroon",
            "mediumaquamarine",
            "mediumblue",
            "mediumorchid",
            "mediumpurple",
            "mediumseagreen",
            "mediumslateblue",
            "mediumspringgreen",
            "mediumturquoise",
            "mediumvioletred",
            "midnightblue",
            "mintcream",
            "mistyrose",
            "moccasin",
            "navajowhite",
            "navy",
            "oldlace",
            "olive",
            "olivedrab",
            "orange",
            "orangered",
            "orchid",
            "palegoldenrod",
            "palegreen",
            "paleturquoise",
            "palevioletred",
            "papayawhip",
            "peachpuff",
            "peru",
            "pink",
            "plum",
            "powderblue",
            "purple",
            "red",
            "rosybrown",
            "royalblue",
            "saddlebrown",
            "salmon",
            "sandybrown",
            "seagreen",
            "seashell",
            "sienna",
            "silver",
            "skyblue",
            "slateblue",
            "slategrey",
            "snow",
            "springgreen",
            "steelblue",
            "tan",
            "teal",
            "thistle",
            "tomato",
            "turquoise",
            "violet",
            "wheat",
            "white",
            "whitesmoke",
            "yellow",
            "yellowgreen"
    };
    std::cout<<std::endl<<std::endl<<std::endl<<"============GFA DUMP STARTING============"<<std::endl;
    std::cout<<"Graph has "<< hb.EdgeObjectCount() <<" edges"<<std::endl;

    // define the k-1 overlap
    int overlap_length = int(hb.K()) - 1;

    // remove any path from filename
    int pos = filename.find_last_of("/");
    std::string prefix = filename;
    if (pos != std::string::npos){
    	prefix = filename.substr(pos + 1);
    }

    vec<vec<vec<vec<int>>>> lines;
    vec<int> to_left, to_right;
    hb.ToLeft(to_left), hb.ToRight(to_right);
    std::vector<int64_t> colour(hb.EdgeObjectCount(), -1);

    if (find_lines) {
        Ofstream(gfa_out, filename + "_lines.gfa");
        Ofstream(fasta_out, filename + "_lines.fasta");
        FindLines(hb, inv, lines, MAX_CELL_PATHS, MAX_DEPTH);
        SortLines(lines, hb, inv);

        gfa_out << "H\tVN:Z:1.0" << std::endl;
        std::vector<int64_t> canonical_included(hb.EdgeObjectCount(), -1);

        int64_t current_colour = 1;
        //First step, mark Edges as used if they appear in a line
        for (auto line : lines) {
            std::vector<std::pair<uint64_t, bool>> prev_segment_end_edges;
            for (auto segment : line) {//or cell, or bubble
                std::vector<std::pair<uint64_t, bool>> end_edges;
                for (auto path : segment) {//or unitig?-ish
                    if (path.empty()) {//empty path (i.e., gap!)
                        end_edges = prev_segment_end_edges;//HACK to not disconnect
                    }
                    else {
                        int64_t prev_in_path = -1;
                        bool prev_in_path_fw;
                        for (auto edge: path) {
                            if (canonical_included[edge] == -1) {
                                uint64_t ce = edge;
                                auto eo=hb.EdgeObject(edge);
                                if (eo.getCanonicalForm()==CanonicalForm::REV) {
                                    ce = inv[edge];
                                    eo=hb.EdgeObject(ce);
                                }
                                canonical_included[edge] = ce;
                                canonical_included[inv[edge]] = ce;

                                gfa_out << "S\tedge" << ce << "\t*"
                                << "\tLN:i:" << eo.isize()
                                << "\tCL:Z:" << colour_names[current_colour % colour_names.size()]
                                << "\tUR:Z:" << prefix << "_lines.fasta"
                                << std::endl;

                                //gfa_out << "S\tedge" << ce << "\t" << hb.EdgeObject(ce) << "\tCL:Z:" <<
                                //colour_names[current_colour % colour_names.size()] << std::endl;

                                colour[ce] = current_colour;
                                colour[inv[ce]] = current_colour;

                                // write seq to FASTA
                                fasta_out << ">edge" << ce << std::endl << eo << std::endl;

                                //std::cout<<"edge"<<ce<<": "<<colour_names[current_colour%colour_names.size()]<<std::endl;
                            } //else {
                            //colour[canonical_included[edge]]=0;
                            //std::cout<<"edge"<<canonical_included[edge]<<": black"<<std::endl;
                            //}
                            if (prev_in_path != -1) {
                                gfa_out << "L\tedge" << prev_in_path << (prev_in_path_fw ? "\t+\tedge" : "\t-\tedge") <<
                                canonical_included[edge] << ((canonical_included[edge] == edge) ? "\t+" : "\t-") <<
                                "\t" << overlap_length << "M" << std::endl;
                            }
                            prev_in_path = canonical_included[edge];
                            prev_in_path_fw = (canonical_included[edge] == edge);

                        }
                        //Connect all previous elements to the first element in the path
                        uint64_t ce = canonical_included[path[0]];
                        bool ce_fw = (ce == path[0]);
                        for (auto pe:prev_segment_end_edges) {
                            gfa_out << "L\tedge" << pe.first << (pe.second ? "\t+\tedge" : "\t-\tedge") << ce <<
                            (ce_fw ? "\t+" : "\t-") << "\t" << overlap_length << "M" << std::endl;
                        }
                        //add last element in the path to the end_elements
                        end_edges.push_back(std::make_pair(prev_in_path, prev_in_path_fw));
                    }
                }
                prev_segment_end_edges = end_edges;

            }
            ++current_colour;
        }
    }
    for (auto &e:marked_edges) {
        colour[e]=112;
        if (inv[e]>=0) colour[inv[e]]=112;
    }
    // Do the raw dump
    Ofstream(gfa_raw_out,filename+"_raw.gfa");
    Ofstream(fasta_raw_out,filename+"_raw.fasta");
    std::cout<<"Dumping edges"<<std::endl;

    // header line
    gfa_raw_out << "H\tVN:Z:1.0" << std::endl;

    for (auto ei=0;ei<hb.EdgeObjectCount();++ei){
        auto eo=hb.EdgeObject(ei);
        if (eo.getCanonicalForm()==CanonicalForm::REV) continue;

        gfa_raw_out << "S\tedge" << ei << "\t*"
        << "\tLN:i:" << eo.isize()
        << "\tCL:Z:" << (colour[ei]>0 ? colour_names[colour[ei]%colour_names.size()] : "black" )
        << "\tUR:Z:" << prefix << "_raw.fasta"
        << std::endl;

        // write seq to FASTA
        fasta_raw_out << ">edge" << ei << std::endl << eo << std::endl;
    }

    std::cout<<"Dumping connections"<<std::endl;
    std::vector<std::vector<uint64_t>> next_edges,prev_edges;
    //TODO: this is stupid duplication of the digraph class, but it's so weird!!!
    prev_edges.resize(to_left.size());
    next_edges.resize(to_right.size());
    for (auto e=0;e<to_left.size();++e){

        uint64_t prev_node=to_left[e];

        prev_edges[e].resize(hb.ToSize(prev_node));
        for (int i=0;i<hb.ToSize(prev_node);++i){
            prev_edges[e][i]=hb.EdgeObjectIndexByIndexTo(prev_node,i);
        }

        uint64_t next_node=to_right[e];

        next_edges[e].resize(hb.FromSize(next_node));
        for (int i=0;i<hb.FromSize(next_node);++i){
            next_edges[e][i]=hb.EdgeObjectIndexByIndexFrom(next_node,i);
        }
    }

    for (uint64_t e=0;e<next_edges.size();e++) {
        //only process the canonical edge
        auto eo=hb.EdgeObject(e);
        if (eo.getCanonicalForm()==CanonicalForm::REV) continue;

        std::set<uint64_t> all_next;
        all_next.insert(next_edges[e].begin(),next_edges[e].end());
        for (auto pi:prev_edges[inv[e]]) all_next.insert(inv[pi]);

        for (auto n:all_next){
            //only process if the canonical of the connection is greater (i.e. only processing "canonical connections")

            uint64_t cn=(hb.EdgeObject(n).getCanonicalForm()!=CanonicalForm ::REV ? n:inv[n]);
            if (cn<e) continue;
            gfa_raw_out << "L\tedge" << e << "\t+\tedge" << cn << (cn==n ? "\t+" : "\t-") << "\t" << overlap_length << "M" << std::endl;
        }

        std::set<uint64_t> all_prev;
        all_prev.insert(prev_edges[e].begin(),prev_edges[e].end());
        for (auto ni:next_edges[inv[e]]) all_prev.insert(inv[ni]);

        for (auto p:all_prev){
            //only process if the canonical of the connection is greater (i.e. only processing "canonical connections")
            uint64_t cp=(hb.EdgeObject(p).getCanonicalForm()!=CanonicalForm ::REV?p:inv[p]);
            if (cp<e) continue;
            gfa_raw_out << "L\tedge" << e << "\t-\tedge" << cp << (cp==p ? "\t-" : "\t+") << "\t" << overlap_length << "M" << std::endl;
        }
    }


    std::cout<<"============GFA DUMP ENDED============"<<std::endl<<std::endl<<std::endl<<std::endl;

}


void GFADumpDetail (std::string filename, const HyperBasevector &hb, const vec<int> &inv, const std::vector<int> & marked_edges, const std::vector<int> & marked_vertices){
    std::vector<std::string> colour_names={
            "aliceblue",
            "antiquewhite",
            "aqua",
            "aquamarine",
            "azure",
            "beige",
            "bisque",
            "blanchedalmond",
            "blue",
            "blueviolet",
            "brown",
            "burlywood",
            "cadetblue",
            "chartreuse",
            "chocolate",
            "coral",
            "cornflowerblue",
            "cornsilk",
            "crimson",
            "cyan",
            "darkblue",
            "darkcyan",
            "darkgoldenrod",
            "darkgreen",
            "darkgrey",
            "darkkhaki",
            "darkmagenta",
            "darkolivegreen",
            "darkorange",
            "darkorchid",
            "darkred",
            "darksalmon",
            "darkseagreen",
            "darkslateblue",
            "darkslategrey",
            "darkturquoise",
            "darkviolet",
            "deeppink",
            "deepskyblue",
            "dimgrey",
            "dodgerblue",
            "firebrick",
            "floralwhite",
            "forestgreen",
            "fuchsia",
            "gainsboro",
            "ghostwhite",
            "gold",
            "goldenrod",
            "grey",
            "green",
            "greenyellow",
            "honeydew",
            "hotpink",
            "indianred",
            "indigo",
            "ivory",
            "khaki",
            "lavender",
            "lavenderblush",
            "lawngreen",
            "lemonchiffon",
            "lightblue",
            "lightcoral",
            "lightcyan",
            "lightgoldenrodyellow",
            "lightgreen",
            "lightgrey",
            "lightpink",
            "lightsalmon",
            "lightseagreen",
            "lightskyblue",
            "lightslategrey",
            "lightsteelblue",
            "lightyellow",
            "lime",
            "limegreen",
            "linen",
            "magenta",
            "maroon",
            "mediumaquamarine",
            "mediumblue",
            "mediumorchid",
            "mediumpurple",
            "mediumseagreen",
            "mediumslateblue",
            "mediumspringgreen",
            "mediumturquoise",
            "mediumvioletred",
            "midnightblue",
            "mintcream",
            "mistyrose",
            "moccasin",
            "navajowhite",
            "navy",
            "oldlace",
            "olive",
            "olivedrab",
            "orange",
            "orangered",
            "orchid",
            "palegoldenrod",
            "palegreen",
            "paleturquoise",
            "palevioletred",
            "papayawhip",
            "peachpuff",
            "peru",
            "pink",
            "plum",
            "powderblue",
            "purple",
            "red",
            "rosybrown",
            "royalblue",
            "saddlebrown",
            "salmon",
            "sandybrown",
            "seagreen",
            "seashell",
            "sienna",
            "silver",
            "skyblue",
            "slateblue",
            "slategrey",
            "snow",
            "springgreen",
            "steelblue",
            "tan",
            "teal",
            "thistle",
            "tomato",
            "turquoise",
            "violet",
            "wheat",
            "white",
            "whitesmoke",
            "yellow",
            "yellowgreen"
    };

    // define the k-1 overlap
    int overlap_length = int(hb.K()) - 1;

    // remove any path from filename
    int pos = filename.find_last_of("/");
    std::string prefix = filename;
    if (pos != std::string::npos){
        prefix = filename.substr(pos + 1);
    }

    vec<vec<vec<vec<int>>>> lines;
    vec<int> to_left, to_right;
    hb.ToLeft(to_left), hb.ToRight(to_right);
    std::vector<int64_t> edge_colour(hb.EdgeObjectCount(), -1);
    std::vector<int64_t> vertex_colour(hb.N(), -1);
    for (auto &e:marked_edges) {
        auto eo=hb.EdgeObject(e);
        if (eo.getCanonicalForm()==CanonicalForm::REV) {
            edge_colour[e] = 420;//pink for reverse edges
        } else {
            edge_colour[e] = 112;
        }
    }
    for (auto &e:marked_vertices) {
        vertex_colour[e] = 423;//pink for reverse edges
    }

    // Do the raw dump
    Ofstream(gfa_raw_out,filename+"_raw.gfa");
    Ofstream(fasta_raw_out,filename+"_raw.fasta");

    // header line
    gfa_raw_out << "H\tVN:Z:1.0" << std::endl;

    for (auto ei=0;ei<hb.EdgeObjectCount();++ei){
        auto eo=hb.EdgeObject(ei);

        gfa_raw_out << "S\tedge" << ei << "\t*"
                    << "\tLN:i:" << eo.isize()
                    << "\tCL:Z:" << (edge_colour[ei]>0 ? colour_names[edge_colour[ei]%colour_names.size()] : (eo.getCanonicalForm()==CanonicalForm::REV ? "grey" : "black" ))
                    << "\tUR:Z:" << prefix << "_raw.fasta"
                    << std::endl;

        // write seq to FASTA
        fasta_raw_out << ">edge" << ei << std::endl << eo << std::endl;
    }

    for (auto vi=0;vi<hb.N();++vi) {
        gfa_raw_out << "S\tvertex" << vi << "\t*"
                    << "\tLN:i:0"
                    << "\tCL:Z:" << (vertex_colour[vi]>0 ? "magenta" : "blue" )
                    << std::endl;
    }
    for (auto vi=0;vi<hb.N();++vi) {
        for (auto fi=0;fi<hb.From(vi).isize();++fi){
            auto ei=hb.EdgeObjectIndexByIndexFrom(vi,fi);
            gfa_raw_out << "L\tvertex" << vi << "\t+\tedge" << ei << "\t+\t0M"<< std::endl;
        }
        for (auto ti=0;ti<hb.To(vi).isize();++ti){
            auto ei=hb.EdgeObjectIndexByIndexTo(vi,ti);
            gfa_raw_out << "L\tedge" << ei << "\t+\tvertex" << vi << "\t+\t0M"<< std::endl;
        }
    }

}

bool check_from_to(const HyperBasevector &hb){
    vec<int> toLeft,toRight;
    hb.ToLeft(toLeft);
    hb.ToRight(toRight);
    bool ok=true;
    for (auto vi=0;vi<hb.N();++vi){
        for (auto fi=0;fi<hb.From(vi).isize();++fi){
            auto ei=hb.EdgeObjectIndexByIndexFrom(vi,fi);
            if (toLeft[ei]!=vi){
                std::cout<<"Vertex #"<<vi<<" claims edge #"<<ei<<" is from it, but edge starts at vertex #"<<toLeft[ei]<<std::endl;
                ok=false;
            }
        }
        for (auto ti=0;ti<hb.To(vi).isize();++ti){
            auto ei=hb.EdgeObjectIndexByIndexTo(vi,ti);
            if (toRight[ei]!=vi){
                std::cout<<"Vertex #"<<vi<<" claims edge #"<<ei<<" is To it, but edge ends at vertex #"<<toRight[ei]<<std::endl;
                ok=false;
            }
        }
    }
    return ok;
}