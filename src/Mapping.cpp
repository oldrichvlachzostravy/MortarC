/*
 * Mapping.cpp
 *
 *  Created on: Mar 19, 2013
 *      Author: Oldřich Vlach
 */


#include "Mapping.h"

void SegmentTriangle::write_ostream(std::ostream& stream, int counter) const
{
	stream  << "    Tri:  " << counter << "      s[ ("
			<< std::setw(7) << this->s[0].x << ", "
			<< std::setw(7) << this->s[0].y << "), ("
			<< std::setw(7) << this->s[1].x << ", "
			<< std::setw(7) << this->s[1].y << "), ("
			<< std::setw(7) << this->s[2].x << ", "
			<< std::setw(7) << this->s[2].y << ")],   m[ ("
			<< std::setw(7) << this->m[0].x << ", "
			<< std::setw(7) << this->m[0].y << "), ("
			<< std::setw(7) << this->m[1].x << ", "
			<< std::setw(7) << this->m[1].y << "), ("
			<< std::setw(7) << this->m[2].x << ", "
			<< std::setw(7) << this->m[2].y << ")]\n";
}

void SegmentLine::write_ostream(std::ostream& stream, int counter) const
{
	stream  << "    Lin:  " << counter << "      s[ "
			<< std::setw(7) << this->s[0] << ", "
			<< std::setw(7) << this->s[1] << "]  ---   m[ "
			<< std::setw(7) << this->m[0] << ", "
			<< std::setw(7) << this->m[1] << " ]\n";
}

template <>
void Mappings<SegmentLine>::mortar_segmentation(Element *element, std::map<int, std::vector<Element*> > &adjacentMap)
{
    // TODO
	double line_projection_matrix[6];
	element->compute_line_projection_matrix(line_projection_matrix);
	element->project_element_to_line(line_projection_matrix);
	Element *e; // master element
	ElementMap *incident_elms = new ElementMap();
	// for each slave element node
	for (int i = 0; i < element->get_node_count(); i++) {
		// find master element that belongs to i-th node
		e = element->get_closest_element(i);
		if (e != NULL) {
			e->project_element_to_line(line_projection_matrix);
			incident_elms->insert(std::pair<int, Element*>(e->get_id(), e));
			// recursively enrich this set by neighboring master elements as much as possible (and correct)
			Element::add_incident_elements_in_line(adjacentMap, e, element, incident_elms, line_projection_matrix);
		}
	} // incident_elms now contains all needed master elements

	SegmentLine reference, absolute;
	MCVec3 u, v, w, t, tp;
	PointList *polygon;
	std::vector<MCVec2> given_projected_nodes;
	FEPrimalBase fe(1);
	mappings.push_back(Mapping<SegmentLine> (element));
	for(ElementMap::iterator it = incident_elms->begin(); it != incident_elms->end(); it++)
	{
		polygon = Element::clip_element_polygons(it->second, element);
		if(polygon->size() == 2)
		{
			//polygons->insert(std::pair<int, PointList*>(it->first, polygon));
			e = it->second;    // master element
			absolute  = SegmentLine();
			reference = SegmentLine();
			given_projected_nodes.resize(2);
			given_projected_nodes[0] = MCVec2(polygon->begin()->x,0);
			given_projected_nodes[1] = MCVec2(polygon->rbegin()->x,0);
			const std::vector<MCVec2> result_slave  = fe.get_reference_coordinates(element, &given_projected_nodes);
			const std::vector<MCVec2> result_master = fe.get_reference_coordinates(      e, &given_projected_nodes);
			double reference_element_area;
			switch (element->get_type())
			{
			case M_ELEMENT_LINE2:
			case M_ELEMENT_LINE3: reference_element_area = 2.0; break;
			default:              fprintf(stderr, "Mappings<SegmentLine>::mortar_segmentation ... called with 3D element\n"); exit(1); break;
			}
			double area = fabs(result_slave[1].x-result_slave[0].x)/reference_element_area;
			mappings.rbegin()->add_to_element_coverage_area_ratio(area);
			SegmentLine segment_line;
			segment_line.m[0] = result_master[0].x;  segment_line.m[1] = result_master[1].x;
			segment_line.s[0] =  result_slave[0].x;  segment_line.s[1] =  result_slave[1].x;
//			std::cout << "S:" << element->get_id() << "   M:" << e->get_id() << "   " <<
//					"xiS=[" << segment_line.s[0] << "," << segment_line.s[1] << "]  " <<
//					"xiM=[" << segment_line.m[0] << "," << segment_line.m[1] << "]  " << std::endl;
//			for (int j = 0; j < 2; j++) {
//				if ( (segment_line.s[j] == 1.0)||(segment_line.s[j] == -1.0) ) {
//					// according to (44) change segment_line.m[j]
//					double xi_k = segment_line.m[j];
//					FEPrimalBase fe_master(3);
//					fe_master.init_all(e);
//					std::vector<MCVec2> master_ref_points(1);
//					for (int k = 0; k < 15; k++) {
//						master_ref_points[0] = xi_k;
//						fe_master.init_reference(e, &master_ref_points);
//						double tmpB = -(element->get_node(j)->get_coordinates().x * element->get_node(j)->get_normal().y - element->get_node(j)->get_coordinates().y * element->get_node(j)->get_normal().x);
//						double tmpA = 0.0;
//						for (int l = 0; l < e->get_node_count(); l++) {
//							tmpB +=  fe_master.get_n()[l][0] * //xl_x_nj[l];
//									e->get_node(l)->get_coordinates().x * element->get_node(j)->get_normal().y - e->get_node(l)->get_coordinates().y * element->get_node(j)->get_normal().x;
//							tmpA +=  fe_master.get_dndxi()[l][0][0] * //xl_x_nj[l];
//									e->get_node(l)->get_coordinates().x * element->get_node(j)->get_normal().y - e->get_node(l)->get_coordinates().y * element->get_node(j)->get_normal().x;
//						}
//						xi_k += - tmpB/tmpA;
//					}
//					if (xi_k < -1.0) xi_k = -1.0;
//					if (xi_k >  1.0) xi_k =  1.0;
//					segment_line.m[j] = xi_k;
//				}
//				if ( (segment_line.m[j] == 1.0)||(segment_line.m[j] == -1.0) ) {
//					// according to (45) change segment_line.s[i]
//					double xi_k = segment_line.s[j];
//					FEPrimalBase fe_slave(3);
//					fe_slave.init_all(e);
//					std::vector<MCVec2> slave_ref_points(1);
//					for (int k = 0; k < 15; k++) {
//						slave_ref_points[0] = xi_k;
//						fe_slave.init_reference(element, &slave_ref_points);
//						MCVec2 tmpvA = MCVec2( 0.0, 0.0);
//						MCVec2 tmpvB = MCVec2( 0.0, 0.0);
//						MCVec2 tmpvC = MCVec2( 0.0, 0.0);
//						MCVec2 tmpvD = MCVec2( 0.0, 0.0);
//						MCVec2 x_l = MCVec2(e->get_node(j)->get_coordinates().x,e->get_node(j)->get_coordinates().y);
//						for (int j = 0; j < element->get_node_count(); j++) {
//							tmpvA += MCVec2(element->get_node(j)->get_coordinates().x,element->get_node(j)->get_coordinates().y) * fe_slave.get_dndxi()[j][0][0];
//							tmpvB += MCVec2(element->get_node(j)->get_coordinates().x,element->get_node(j)->get_coordinates().y) * fe_slave.get_n()[j][0];
//							tmpvC += MCVec2(element->get_node(j)->get_normal().x,element->get_node(j)->get_normal().y) * fe_slave.get_dndxi()[j][0][0];
//							tmpvC += MCVec2(element->get_node(j)->get_normal().x,element->get_node(j)->get_normal().y) * fe_slave.get_n()[j][0];
//						}
//						double tmpA = tmpvA.x*tmpvC.y - tmpvA.y*tmpvC.x;
//						double tmpB = -( (tmpvB.x-x_l.x)*tmpvD.y - (tmpvB.y-x_l.y)*tmpvD.x );
//						xi_k += -tmpB/tmpA;
//					}
//					if (xi_k < -1.0) xi_k = -1.0;
//					if (xi_k >  1.0) xi_k =  1.0;
//					segment_line.s[j] = xi_k;
//				}
//
//			}
//			std::cout << "S:" << element->get_id() << "   M:" << e->get_id() << "   " <<
//								"xiS=[" << segment_line.s[0] << "," << segment_line.s[1] << "]  " <<
//								"xiM=[" << segment_line.m[0] << "," << segment_line.m[1] << "]  " << std::endl;
//			std::cout << std::endl;
			mappings.rbegin()->add_segments_for_master(e, segment_line);
		}
		delete polygon;
	}
    delete incident_elms;
}

template <>
void Mappings<SegmentTriangle>::mortar_segmentation(Element *element, std::map<int, std::vector<Element*> > &adjacentMap)
{
    double plane_projection_matrix[12];
    element->compute_plane_projection_matrix(plane_projection_matrix);
    element->project_element_to_plane(plane_projection_matrix);


    Element *e;
    ElementMap *incident_elms = new ElementMap();
    for (int i = 0; i < element->get_node_count(); i++) {
        e = element->get_closest_element(i);
        if (e != NULL) {
            e->project_element_to_plane(plane_projection_matrix);
            // for each slave element node search for master elements belonging to
            incident_elms->insert(std::pair<int, Element*>(e->get_id(), e));
            // recursively enrich this set by neighboring master elements as much as possible (and correct)
            Element::add_incident_elements_in_plane(adjacentMap, e, element, incident_elms, plane_projection_matrix);
        }
    } // incident_elms now contains all needed master elements

    PolygonMap *polygons = new PolygonMap();
    PointList *polygon;
    for(ElementMap::iterator it = incident_elms->begin(); it != incident_elms->end(); it++)
    {
        polygon = Element::clip_element_polygons(it->second, element);
        if(polygon->size() > 2)
        {
            polygons->insert(std::pair<int, PointList*>(it->first, polygon));
        }
        else
        {
            delete polygon;
        }
    }
    SegmentTriangle reference, absolute;
    MCVec3 u, v, w, t, tp;

    std::vector<TDescription*> triangles;
    PolygonMap::iterator pit;
    mappings.push_back(Mapping<SegmentTriangle> (element));
    // decompose all polygons to triangles
    for (pit = polygons->begin(); pit != polygons->end(); pit++) {
        e = incident_elms->at(pit->first);    // master element
        triangles = triangulate(pit->second);
        FEPrimalBase fe(1);
        std::vector<MCVec2> given_projected_nodes;
        for (uint j = 0; j < triangles.size(); j++) {
            absolute = SegmentTriangle();
            reference = SegmentTriangle();
            given_projected_nodes.resize(3);
            given_projected_nodes[0] = triangles[j]->p[0];
            given_projected_nodes[1] = triangles[j]->p[1];
            given_projected_nodes[2] = triangles[j]->p[2];
            const std::vector<MCVec2> result_slave  = fe.get_reference_coordinates(element, &given_projected_nodes);
            const std::vector<MCVec2> result_master = fe.get_reference_coordinates(e, &given_projected_nodes);

            // area of triangle ( see http://www.mathopenref.com/coordtrianglearea.html)
            double reference_element_area;
            switch (element->get_type())
            {
            case M_ELEMENT_LINE2:
            case M_ELEMENT_LINE3: reference_element_area = 2.0; break;
            case M_ELEMENT_TRIA3:
            case M_ELEMENT_TRIA6: reference_element_area = 0.5; break;
            case M_ELEMENT_QUAD4:
            case M_ELEMENT_QUAD8: reference_element_area = 4.0; break;
            default:              reference_element_area = 1.0; break;
            }
            double area = 0.5*fabs(
                    result_slave[0].x*(result_slave[1].y-result_slave[2].y) +
                    result_slave[1].x*(result_slave[2].y-result_slave[0].y) +
                    result_slave[2].x*(result_slave[0].y-result_slave[1].y))/reference_element_area;
            mappings.rbegin()->add_to_element_coverage_area_ratio(area);
            SegmentTriangle segment_triangle;
            segment_triangle.m[0] = result_master[0];  segment_triangle.m[1] = result_master[1];  segment_triangle.m[2] = result_master[2];
            segment_triangle.s[0] =  result_slave[0];  segment_triangle.s[1] =  result_slave[1];  segment_triangle.s[2] =  result_slave[2];
            mappings.rbegin()->add_segments_for_master(e, segment_triangle);
            delete triangles[j];
        }
    }

    delete incident_elms;
    for(pit = polygons->begin(); pit != polygons->end(); pit++) {
        delete pit->second;
    }
    delete polygons;
}

template <>
int Mapping<SegmentLine>::write_ensight_gold( std::ofstream *geofile_ptr, int &node_counter, int &element_counter)
{
    if (geofile_ptr->is_open())
    {
        int tmp_int = 0;
        int tmp_int_counter=1;
        /* COORDINATES */
        (*geofile_ptr) << "coordinates\n";
        for (std::map<int, std::vector<SegmentLine> >::iterator sit = segments_for_master.begin(); sit != segments_for_master.end(); ++sit)
        {
            tmp_int += sit->second.size();
        }
        (*geofile_ptr) << std::setw( ENSIGHT_GOLD_INT_WIDTH) << 2*tmp_int << "\n";
        FEPrimalBase fe(1);
        for (int i=0; i<2*tmp_int; i++)
        {
            (*geofile_ptr) << std::setw( ENSIGHT_GOLD_INT_WIDTH) << node_counter++ << "\n";
        }
        std::vector<MCVec3> tmp_vector_mcvec3;
        for (std::map<int, std::vector<SegmentLine> >::iterator sit= segments_for_master.begin(); sit != segments_for_master.end(); ++sit)
        {
            for (unsigned int i=0; i < sit->second.size(); i++)
            {
                std::vector<MCVec2> tmp_vector_mcvec2;
                tmp_vector_mcvec2.push_back(MCVec2(sit->second[i].s[0]));
                tmp_vector_mcvec2.push_back(MCVec2(sit->second[i].s[1]));
                fe.init_all(element_slave,&tmp_vector_mcvec2);
                const std::vector<MCVec3> tmp_mcvec3 = fe.get_refpoints_coordiantes(&tmp_vector_mcvec2);
                tmp_vector_mcvec3.push_back(tmp_mcvec3[0]);
                tmp_vector_mcvec3.push_back(tmp_mcvec3[1]);
            }
        }
        for (std::vector<MCVec3>::iterator it = tmp_vector_mcvec3.begin(); it != tmp_vector_mcvec3.end(); ++it)
        {
            (*geofile_ptr) << std::setw( ENSIGHT_GOLD_DOUBLE_WIDTH) << it->x << "\n";
        }
        for (std::vector<MCVec3>::iterator it = tmp_vector_mcvec3.begin(); it != tmp_vector_mcvec3.end(); ++it)
        {
            (*geofile_ptr) << std::setw( ENSIGHT_GOLD_DOUBLE_WIDTH) << it->y << "\n";
        }
        for (std::vector<MCVec3>::iterator it = tmp_vector_mcvec3.begin(); it != tmp_vector_mcvec3.end(); ++it)
        {
            (*geofile_ptr) << std::setw( ENSIGHT_GOLD_DOUBLE_WIDTH) << it->z << "\n";
        }
        /* ELEMENTS */
        (*geofile_ptr) << "bar2\n";
        (*geofile_ptr) << std::setw( ENSIGHT_GOLD_INT_WIDTH) << tmp_int << "\n";
        for (int i=0; i<tmp_int; i++)
        {
            (*geofile_ptr) << std::setw( ENSIGHT_GOLD_INT_WIDTH) << element_counter++<< "\n";
        }
        for (int i=0; i<tmp_int; i++)
        {
            (*geofile_ptr)
                    << std::setw( ENSIGHT_GOLD_INT_WIDTH) << tmp_int_counter
                    << std::setw( ENSIGHT_GOLD_INT_WIDTH) << tmp_int_counter+1 << "\n";
            tmp_int_counter += 2;
        }
        return 1;
    }
    else
        return 0;
}

template <>
int Mapping<SegmentTriangle>::write_ensight_gold( std::ofstream *geofile_ptr, int &node_counter, int &element_counter)
{
    if (geofile_ptr->is_open())
    {
        int tmp_int = 0;
        int tmp_int_counter=1;
        /* COORDINATES */
        (*geofile_ptr) << "coordinates\n";
        for (std::map<int, std::vector<SegmentTriangle> >::iterator sit = segments_for_master.begin(); sit != segments_for_master.end(); ++sit)
        {
            tmp_int += sit->second.size();
        }
        (*geofile_ptr) << std::setw( ENSIGHT_GOLD_INT_WIDTH) << 3*tmp_int << "\n";
        FEPrimalBase fe(1);
        for (int i=0; i<3*tmp_int; i++)
        {
            (*geofile_ptr) << std::setw( ENSIGHT_GOLD_INT_WIDTH) << node_counter++ << "\n";
        }
        std::vector<MCVec3> tmp_vector_mcvec3;
        for (std::map<int, std::vector<SegmentTriangle> >::iterator sit= segments_for_master.begin(); sit != segments_for_master.end(); ++sit)
        {
            for (unsigned int i=0; i < sit->second.size(); i++)
            {
                std::vector<MCVec2> tmp_vector_mcvec2;
                tmp_vector_mcvec2.push_back(MCVec2(sit->second[i].s[0]));
                tmp_vector_mcvec2.push_back(MCVec2(sit->second[i].s[1]));
                tmp_vector_mcvec2.push_back(MCVec2(sit->second[i].s[2]));
                fe.init_all(element_slave,&tmp_vector_mcvec2);
                const std::vector<MCVec3> tmp_mcvec3 = fe.get_refpoints_coordiantes(&tmp_vector_mcvec2);
                tmp_vector_mcvec3.push_back(tmp_mcvec3[0]);
                tmp_vector_mcvec3.push_back(tmp_mcvec3[1]);
                tmp_vector_mcvec3.push_back(tmp_mcvec3[2]);
            }
        }
        for (std::vector<MCVec3>::iterator it = tmp_vector_mcvec3.begin(); it != tmp_vector_mcvec3.end(); ++it)
        {
            (*geofile_ptr) << std::setw( ENSIGHT_GOLD_DOUBLE_WIDTH) << it->x << "\n";
        }
        for (std::vector<MCVec3>::iterator it = tmp_vector_mcvec3.begin(); it != tmp_vector_mcvec3.end(); ++it)
        {
            (*geofile_ptr) << std::setw( ENSIGHT_GOLD_DOUBLE_WIDTH) << it->y << "\n";
        }
        for (std::vector<MCVec3>::iterator it = tmp_vector_mcvec3.begin(); it != tmp_vector_mcvec3.end(); ++it)
        {
            (*geofile_ptr) << std::setw( ENSIGHT_GOLD_DOUBLE_WIDTH) << it->z << "\n";
        }
        /* ELEMENTS */
        (*geofile_ptr) << "tria3\n";
        (*geofile_ptr) << std::setw( ENSIGHT_GOLD_INT_WIDTH) << tmp_int << "\n";
        for (int i=0; i<tmp_int; i++)
        {
            (*geofile_ptr) << std::setw( ENSIGHT_GOLD_INT_WIDTH) << element_counter++<< "\n";
        }
        for (int i=0; i<tmp_int; i++)
        {
            (*geofile_ptr)
                    << std::setw( ENSIGHT_GOLD_INT_WIDTH) << tmp_int_counter
                    << std::setw( ENSIGHT_GOLD_INT_WIDTH) << tmp_int_counter+1
                    << std::setw( ENSIGHT_GOLD_INT_WIDTH) << tmp_int_counter+2 << "\n";
            tmp_int_counter += 3;
        }
        return 1;
    }
    else
        return 0;
}
