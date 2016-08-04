//
// Created by Bernardo Clavijo (TGAC) on 07/07/2016.
//

#include <paths/long/large/Lines.h>
#include "GFADump.h"

void GFADump (std::string filename, const HyperBasevector &hb, const vec<int> &inv, const
ReadPathVec &paths, const int MAX_CELL_PATHS, const int MAX_DEPTH, bool find_lines){

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
    vec<vec<vec<vec<int>>>> lines;
    vec<int> to_left, to_right;
    hb.ToLeft(to_left), hb.ToRight(to_right);
    std::vector<int64_t> colour(hb.EdgeObjectCount(), -1);

    if (find_lines) {
        Ofstream(gfa_out, filename + "_lines.gfa");
        FindLines(hb, inv, lines, MAX_CELL_PATHS, MAX_DEPTH);
        SortLines(lines, hb, inv);


        gfa_out << "H\tVN:Z:1.0" << std::endl;
        std::vector<int64_t> canonical_included(hb.EdgeObjectCount(), -1);

        int64_t current_colour = 1;
        //TODO: Dump the overlaps correctly
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
                                if (hb.EdgeObject(inv[edge]) < hb.EdgeObject(edge)) {
                                    ce = inv[edge];
                                }
                                canonical_included[edge] = ce;
                                canonical_included[inv[edge]] = ce;
                                gfa_out << "S\tedge" << ce << "\t" << hb.EdgeObject(ce) << "\tCL:z:" <<
                                colour_names[current_colour % colour_names.size()] << std::endl;
                                colour[ce] = current_colour;
                                colour[inv[ce]] = current_colour;
                                //std::cout<<"edge"<<ce<<": "<<colour_names[current_colour%colour_names.size()]<<std::endl;
                            } //else {
                            //colour[canonical_included[edge]]=0;
                            //std::cout<<"edge"<<canonical_included[edge]<<": black"<<std::endl;
                            //}
                            if (prev_in_path != -1) {
                                gfa_out << "L\tedge" << prev_in_path << (prev_in_path_fw ? "\t+\tedge" : "\t-\tedge") <<
                                canonical_included[edge] << ((canonical_included[edge] == edge) ? "\t+" : "\t-") <<
                                "\t0M" << std::endl;
                            }
                            prev_in_path = canonical_included[edge];
                            prev_in_path_fw = (canonical_included[edge] == edge);

                        }
                        //Connect all previous elements to the first element in the path
                        uint64_t ce = canonical_included[path[0]];
                        bool ce_fw = (ce == path[0]);
                        for (auto pe:prev_segment_end_edges) {
                            gfa_out << "L\tedge" << pe.first << (pe.second ? "\t+\tedge" : "\t-\tedge") << ce <<
                            (ce_fw ? "\t+" : "\t-") << "\t0M" << std::endl;
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
    //Now reconstruct the to-left, to-right, RAW thing
    /*
    struct s_link_edge_to_vertex{
        uint64_t vert;
        uint64_t edge;
        bool in;
    };
    std::set<struct s_link> links;
    for (uint64_t i=0;i<to_left.size();++i){

        struct s_link_edge_to_vertex l;
        if (inv[i]<i) continue;//only process canonical edges
        l.edge;
        links.insert(l);

    };
    for (uint64_t i=0;i<to_right.size();++i){
        struct s_link l;
        l.e1=i;
        l.o1=true;
        if (inv[l.e1]<l.e1) {
            l.e1=inv[l.e1];
            l.o1=false;
        }
        l.e2=to_right[i];
        l.o2=true;
        if (inv[l.e2]<l.e2) {
            l.e2=inv[l.e2];
            l.o2=false;
        }
        //TODO: dist
        l.dist=0;
        links.insert(l);
    };
    for (auto l:links){*/

    Ofstream(gfa_raw_out,filename+"_raw.gfa");
    std::cout<<"Dumping edges"<<std::endl;
    for (auto ei=0;ei<hb.EdgeObjectCount();++ei){
        if (ei>inv[ei]) continue;
        gfa_raw_out << "S\tedge" << ei <<"\t"<< hb.EdgeObject(ei)
        << "\tCL:z:"<< (colour[ei]>0 ? colour_names[colour[ei]%colour_names.size()] : "black" )
         << std::endl;

    }
    std::cout<<"Dumping connections"<<std::endl;
    for (int64_t i=0;i<to_right.size();++i) {
        for (int64_t j=0;j<to_left.size();++j) {
            if (to_left[j]==to_right[i]) {
                int64_t ce1= (i<inv[i]?i:inv[i]);
                int64_t ce2= (j<inv[j]?j:inv[j]);
                gfa_raw_out << "L\tedge" << ce1 << (ce1==i ? "\t+\tedge" : "\t-\tedge") << ce2 <<
                (ce2==j ? "\t+" : "\t-") << "\t0M" << std::endl;
            }
        }
    }


    std::cout<<"============GFA DUMP ENDED============"<<std::endl<<std::endl<<std::endl<<std::endl;

}