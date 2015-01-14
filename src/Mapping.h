
#ifndef MAPPING_H_
#define MAPPING_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <sstream>

#include "SystemIncludes.h"
#include "MCVector.h"
#include "Element.h"
#include "FEPrimalBase.h"
#include "Utils.h"
#include "Boundary.h"
#include "BoundaryMapper.h"

#ifndef ENSIGHT_GOLD_DOUBLE_WIDTH
#define ENSIGHT_GOLD_DOUBLE_WIDTH 12
#endif
#ifndef ENSIGHT_GOLD_INT_WIDTH
#define ENSIGHT_GOLD_INT_WIDTH 10
#endif

class SegmentTriangle
{
	public:
	    MCVec2 m[3];
	    MCVec2 s[3];
	    void write_ostream(std::ostream& stream, int counter) const;
};

class SegmentLine
{
    public:
	    double m[2];
	    double s[2];
	    void write_ostream(std::ostream& stream, int counter) const;
};

class Element;

/// mapping for one slave element
template <class T>
class Mapping
{
	public:
		Mapping(Element *el_slave) {
		    element_slave = el_slave;
            element_coverage_area_ratio = 0.0;;
		}
		//void add_segments_for_master(Element *, const std::vector<MCVec2> &, const std::vector<MCVec2> &);
		void add_segments_for_master(Element *, const T &);

		Element * get_element_slave()  { return element_slave; };
		double get_element_coverage_area_ratio() {return element_coverage_area_ratio; };
		void add_to_element_coverage_area_ratio(double area) {element_coverage_area_ratio += area; };
		const std::map<int, std::vector<T> > & get_segments_for_master() { return segments_for_master; };

		int write_ensight_gold( std::ofstream*, int&, int&);

	private:
		Element * element_slave;
		double element_coverage_area_ratio;
		typename std::map<int, std::vector<T> > segments_for_master;
};

/// vector of mappings ... i.e. for all slave elements
template <class T>
class Mappings
{
    public:
        void compute_mapping(Boundary *);
        std::vector<Mapping<T> > & get_mappings() { return mappings; }
        //const std::vector<Mapping<T> > & get_const_mappings() const { return mappings; }

        void write_mapping(Boundary * master, const char *fname, int ii = 0);

        //int write_ensight_gold(std::ofstream *geofile_ptr, int &node_counter, int &element_counter);

        void write_ensight_gold_slave_master_mapping(BoundaryMapper &project, const char *fname, int ii = 0);

        void write_ensight_gold_normals(BoundaryMapper &project, const char *fname, int ii = 0);

    private:
        void mortar_segmentation(Element *, std::map<int, std::vector<Element*> > &);

        std::vector<Mapping<T> > mappings;
};

template <class T> void Mappings<T>::compute_mapping(Boundary *slave) {
    std::vector<Element*>::iterator it;
    for (it = slave->get_elements().begin(); it != slave->get_elements().end(); it++) {
        mortar_segmentation(*it, slave->get_adjacent());
    }
}


//template <class T> void Mapping<T>::add_segments_for_master(Element * master_el, const std::vector<MCVec2> & slave, const std::vector<MCVec2> & master)
//{
//    SegmentTriangle tri;
//    tri.m[0] = master[0];  tri.m[1] = master[1];  tri.m[2] = master[2];
//    tri.s[0] =  slave[0];  tri.s[1] =  slave[1];  tri.s[2] =  slave[2];
//    if (segments_for_master.count(master_el->get_id()) < 1)
//    {
//        segments_for_master.insert(std::pair<int, std::vector<SegmentTriangle> >(master_el->get_id(), std::vector<SegmentTriangle>()));
//    }
//    segments_for_master[master_el->get_id()].push_back(tri);
//}
template <class T> void Mapping<T>::add_segments_for_master(Element * master_el, const T & segment)
{
    if (segments_for_master.count(master_el->get_id()) < 1)
    {
        segments_for_master.insert(std::pair<int, std::vector<T> >(master_el->get_id(), std::vector<T>()));
    }
    segments_for_master[master_el->get_id()].push_back(segment);
}



template <class T>
void Mappings<T>::write_mapping(Boundary * master, const char *fname, int ii)
{
    std::ostringstream tmp_ostringstream;
    time_t rawtime;
    struct tm * timeinfo;
    tmp_ostringstream << fname << "/mapping_" << ii << ".txt";
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    std::ofstream out_file;
    out_file.open(tmp_ostringstream.str().c_str());
    if (out_file.is_open())
    {
        out_file << "Mapping description (" << timeinfo->tm_mday << "-" << 1+timeinfo->tm_mon << "-" << 1900+timeinfo->tm_year << ")\n\n";
        out_file << std::setiosflags(std::ios::right) << std::setprecision(5);
        for (size_t i = 0; i < mappings.size(); i++)
        {
            out_file << "S:  " << mappings[i].get_element_slave()->get_id() << " (cr = " << mappings[i].get_element_coverage_area_ratio() << ")   [\n";
            for (int j = 0; j < mappings[i].get_element_slave()->get_node_count(); j++)
            {
                out_file << "      "
                		<< std::setw(4) << mappings[i].get_element_slave()->get_node(j)->get_id() << "  ("
                		<< std::setw( ENSIGHT_GOLD_DOUBLE_WIDTH) << mappings[i].get_element_slave()->get_node(j)->get_coordinates().x << ",  "
                        << std::setw( ENSIGHT_GOLD_DOUBLE_WIDTH) << mappings[i].get_element_slave()->get_node(j)->get_coordinates().y << ",  "
                        << std::setw( ENSIGHT_GOLD_DOUBLE_WIDTH) << mappings[i].get_element_slave()->get_node(j)->get_coordinates().z << " )\n";
            }
            out_file << "               ] -----------------------------------\n";
            for (typename std::map<int, std::vector<T> >::const_iterator it = mappings[i].get_segments_for_master().begin(); it != mappings[i].get_segments_for_master().end(); ++it)
            {
                out_file << "  M:  " << it->first << "  [\n";
                for (int j = 0; j < master->get_element(it->first)->get_node_count(); j++)
                {
                    out_file << "        "
                    		<< std::setw(4) << master->get_element(it->first)->get_node(j)->get_id() << "  ("
                    		<< std::setw( ENSIGHT_GOLD_DOUBLE_WIDTH) << master->get_element(it->first)->get_node(j)->get_coordinates().x << ",  "
                            << std::setw( ENSIGHT_GOLD_DOUBLE_WIDTH) << master->get_element(it->first)->get_node(j)->get_coordinates().y << ",  "
                            << std::setw( ENSIGHT_GOLD_DOUBLE_WIDTH) << master->get_element(it->first)->get_node(j)->get_coordinates().z << " )\n";
                }
                out_file << " ]\n";
                int counter=0;
                for (typename std::vector<T>::const_iterator seg_it = it->second.begin(); seg_it != it->second.end(); ++seg_it)
                {
                    seg_it->write_ostream(out_file, ++counter);
                }
            }
        }
    }
}

template <class T> void Mappings<T>::write_ensight_gold_slave_master_mapping(BoundaryMapper &project, const char *fname, int ii)
{
    std::ostringstream tmp_ostringstream;
    time_t rawtime;
    struct tm * timeinfo;
    int part_counter = 1;
    int node_counter = 1;
    int element_counter = 1;
    int tmp_int;
    // ENSIGHT .case FILE
    tmp_ostringstream << fname << "/output_mapping_" << ii << ".case";
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    std::ofstream out_file;
    out_file.open(tmp_ostringstream.str().c_str());
    std::cout << "filename: " << tmp_ostringstream.str().c_str() << std::endl;
    if (out_file.is_open())
    {
        out_file << "# Date (" << timeinfo->tm_mday << "-" << 1+timeinfo->tm_mon << "-" << 1900+timeinfo->tm_year << ")\n";
        out_file << "# EnSight Gold Model\n\n";
        out_file << "FORMAT\n";
        out_file << "type:                   ensight gold\n\n";
        out_file << "GEOMETRY\n";
        out_file << "model:                  output_mapping_" << "000" << ".geo\n\n";
        out_file.close();
    }
    // ENSIGHT .geo FILE
    tmp_ostringstream.str("");
    tmp_ostringstream << fname << "/output_mapping_" << "000" << ".geo";
    out_file.open(tmp_ostringstream.str().c_str());
    if (out_file.is_open())
    {
        out_file << "This is the description of the EnSight Gold geometry, MODEL:\n";
        out_file << "Date (" << timeinfo->tm_mday << "-" << 1+timeinfo->tm_mon << "-" << 1900+timeinfo->tm_year << ")\n";
        out_file << "node id given\n";
        out_file << "element id given\n";
        out_file << "extents\n";
        out_file << std::setiosflags(std::ios::right | std::ios::scientific) << std::setprecision(5);
        BoundingVolume * tmp_interval = project.get_master_bvt()->get_volume();
        out_file << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(0).start
                 << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(0).end << "\n";
        out_file << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(1).start
                 << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(1).end << "\n";
        out_file << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(4).start
                 << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(4).end << "\n";
        // SLAVE
        out_file << "part\n";
        out_file << std::setw(ENSIGHT_GOLD_INT_WIDTH) << part_counter++ << "\n";
        out_file << "slave\n";
        project.get_slave()->write_ensight_gold(&out_file, node_counter, element_counter);
        // MASTER
        out_file << "part\n";
        out_file << std::setw( ENSIGHT_GOLD_INT_WIDTH) << part_counter++ << "\n";
        out_file << "master\n";
        project.get_master()->write_ensight_gold(&out_file, node_counter, element_counter);
        // MAPPING
        for (size_t i = 0; i< mappings.size(); i++)
        {
            if (mappings[i].get_element_coverage_area_ratio() < MIN_SLAVE_COVER_RATIO) continue;
            out_file << "part\n";
            out_file << std::setw( ENSIGHT_GOLD_INT_WIDTH) << part_counter++ << "\n";
            out_file << "mapping[" << i << "] (" << mappings[i].get_element_slave()->get_id() << ")\n";
            mappings[i].write_ensight_gold(&out_file, node_counter, element_counter);
        }
    }
    out_file.close();
}

template <class T> void Mappings<T>::write_ensight_gold_normals(BoundaryMapper &project, const char *fname, int ii)
{
    std::ostringstream tmp_ostringstream;
    time_t rawtime;
    struct tm * timeinfo;
    int part_counter = 1;
    int node_counter = 1;
    int element_counter = 1;
    // ENSIGHT .case FILE
    std::ofstream out_file;
    tmp_ostringstream << fname << "/output_slave_master_" << ii << ".case";
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    out_file.open(tmp_ostringstream.str().c_str());
    if (out_file.is_open())
    {
        out_file << "# Date (" << timeinfo->tm_mday << "-" << 1+timeinfo->tm_mon << "-" << 1900+timeinfo->tm_year << ")\n";
        out_file << "# EnSight Gold Model\n\n";
        out_file << "FORMAT\n";
        out_file << "type:                   ensight gold\n\n";
        out_file << "GEOMETRY\n";
        out_file << "model:                  output_slave_master_" << ii << ".geo\n\n";
        out_file << "VARIABLES\n";
        out_file << "vector per node: 1 1 Normals output_normals_" << ii << ".Nvec\n";
        out_file << "scalar per node: 1 1 Supports output_supports_" << ii << ".Nsca\n";
        out_file.close();
    }
    // ENSIGHT .geo FILE
    tmp_ostringstream.str("");
    tmp_ostringstream << fname << "/output_slave_master_" << ii << ".geo";
    out_file.open(tmp_ostringstream.str().c_str());
    if (out_file.is_open())
    {
        out_file << "This is the description of the EnSight Gold geometry, MODEL:\n";
        out_file << "Date (" << timeinfo->tm_mday << "-" << 1+timeinfo->tm_mon << "-" << 1900+timeinfo->tm_year << ")\n";
        out_file << "node id given\n";
        out_file << "element id given\n";
        out_file << "extents\n";
        out_file << std::setiosflags(std::ios::right | std::ios::scientific) << std::setprecision(5);
        BoundingVolume * tmp_interval = project.get_master_bvt()->get_volume();
        out_file << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(0).start
                 << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(0).end << "\n";
        out_file << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(1).start
                 << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(1).end << "\n";
        out_file << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(4).start
                 << std::setw(ENSIGHT_GOLD_DOUBLE_WIDTH) << tmp_interval->get_bound(4).end << "\n";
        // SLAVE
        out_file << "part\n";
        out_file << std::setw(ENSIGHT_GOLD_INT_WIDTH) << part_counter++ << "\n";
        out_file << "slave\n";
        project.get_slave()->write_ensight_gold(&out_file, node_counter, element_counter);
        // MASTER
        out_file << "part\n";
        out_file << std::setw( ENSIGHT_GOLD_INT_WIDTH) << part_counter++ << "\n";
        out_file << "master\n";
        project.get_master()->write_ensight_gold(&out_file, node_counter, element_counter);
    }
    out_file.close();
    // ENSIGHT contact_3d_mex_withoutclasses__normals.Nvec FILE
    tmp_ostringstream.str("");
    part_counter = 1;
    tmp_ostringstream << fname << "/output_normals_" << ii << ".Nvec";
    out_file.open(tmp_ostringstream.str().c_str());
    if (out_file.is_open())
    {
        out_file << "This is the description of the EnSight Gold geometry, MODEL:\n";
        // SLAVE
        out_file << "part\n";
        out_file << std::setw(ENSIGHT_GOLD_INT_WIDTH) << part_counter++ << "\n";
        project.get_slave()->write_normals_ensight_gold(&out_file);
        // MASTER
        out_file << "part\n";
        out_file << std::setw(ENSIGHT_GOLD_INT_WIDTH) << part_counter++ << "\n";
        project.get_master()->write_normals_ensight_gold(&out_file);
    }
    out_file.close();
    // ENSIGHT contact_3d_mex_withoutclasses__supports.Nsca FILE
    tmp_ostringstream.str("");
    part_counter = 1;
    tmp_ostringstream << fname << "/output_supports_" << ii << ".Nsca";
    out_file.open(tmp_ostringstream.str().c_str());
    if (out_file.is_open())
    {
        out_file << "This is the description of the EnSight Gold geometry, MODEL:\n";
        // SLAVE
        out_file << "part\n";
        out_file << std::setw(ENSIGHT_GOLD_INT_WIDTH) << part_counter++ << "\n";
        project.get_slave()->write_supports_ensight_gold(&out_file);
        // MASTER
        out_file << "part\n";
        out_file << std::setw(ENSIGHT_GOLD_INT_WIDTH) << part_counter++ << "\n";
        project.get_master()->write_supports_ensight_gold(&out_file);
    }
    out_file.close();
}

#endif /* MAPPING_H_ */
