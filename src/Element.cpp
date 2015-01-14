#include "Element.h"
#include "Node.h"
#include "GaussianQuadrature.h"

void Element::compute_center()
{
	MCVec3 c = MCVec3(0, 0, 0);
	for(int i = 0; i < this->node_count; i++) {
		c += this->nodes[i]->get_coordinates();
	}
	c /= this->node_count;
	this->center = new Node(-1, c);
}

void Element::set_bound_volume(BoundingVolume *bounds)
{
	this->element_bounds = bounds;
}

BoundingVolume * Element::get_element_bounds()
{
	return this->element_bounds;
}

Interval Element::get_element_bound(int bound_index)
{
	return this->element_bounds->get_bound(bound_index);
}

void Element::calculate_centers_normal()
{
	MCVec3 n;
	MCVec3 *j = this->get_jacobian(0, 0);
	switch (this->get_type())
	{
	case M_ELEMENT_LINE2:
	case M_ELEMENT_LINE3:
		n = MCVec2(-j->y,j->x); // rotate 90
		n.normalize();
		this->center->set_normal(n);
		break;
	case M_ELEMENT_TRIA3:
	case M_ELEMENT_TRIA6:
	case M_ELEMENT_QUAD4:
	case M_ELEMENT_QUAD8:
		n = cross_prod(j[0], j[1]);
		n.normalize();
		delete[] j;
		this->center->set_normal(n);
		break;
	default:
	    break;
	}
}

int Element::get_id()
{
	return this->id;
}

Node * Element::get_node(int i)
{
	return this->nodes[i];
}

Node ** Element::get_nodes()
{
	return this->nodes;
}

int Element::get_node_count()
{
	return this->node_count;
}

void Element::set_divide_flag(bool flag)
{
	this->divide_flag = flag;
}

bool Element::get_divide_flag()
{
	return this->divide_flag;
}

Element* Element::get_closest_element(int index)
{
    return closest_elements[index];
}

void Element::set_closest_element(Element *e, int index)
{
	this->closest_elements[index] = e;
}

void Element::compute_plane_projection_matrix(double *matrix)
{
	/**
	 *  See E. Sojka: Computer Graphics I - "skripta" http://mrl.cs.vsb.cz/people/sojka/pg/pocitacova_grafikaII.pdf
	 *  p. 22
	 */
	MCVec3 z = center->get_normal();
	MCVec3 x = MCVec3(1, 0, 0);
	MCVec3 y = MCVec3(0, 1, 0);
	MCVec3 t = center->get_coordinates();
	t.flip();
	if(z.x != 0 || z.y != 0) {
		x = MCVec3(-z.y, z.x, 0);
		x.normalize();
		y = cross_prod(z, x);
		y.normalize();
	}
	if(z.z == -1) {
		x = MCVec3(1, 0, 0);
		y = cross_prod(z, x);
		y.normalize();
	}
	matrix[0] = x.x;
	matrix[1] = y.x;
	matrix[2] = z.x;
	matrix[3] = x.y;
	matrix[4] = y.y;
	matrix[5] = z.y;
	matrix[6] = x.z;
	matrix[7] = y.z;
	matrix[8] = z.z;
	matrix[9] = dot_prod(x, t);
	matrix[10] = dot_prod(y, t);
	matrix[11] = dot_prod(z, t);
}

void Element::compute_line_projection_matrix(double *matrix)
{
	MCVec2 y = MCVec2(center->get_normal().x,center->get_normal().y);
    MCVec2 x = MCVec2(1, 0);
	MCVec2 t = MCVec2(center->get_coordinates().x, center->get_coordinates().y);
    t.flip();
	if(y.x != 0) {
		x = MCVec2(-y.y, y.x);
		x.normalize();
	}
	matrix[0] = x.x;
	matrix[1] = y.x;
	matrix[2] = x.y;
	matrix[3] = y.y;
	matrix[4] = dot_prod(x, t);
	matrix[5] = dot_prod(y, t);
}

void Element::project_element_to_plane(double *plane_projection_matrix)
{
	for(int i = 0; i < node_count; i++) {
		nodes[i]->project_point_to_plane(plane_projection_matrix);
	}
}

void Element::project_element_to_line(double *line_projection_matrix)
{
	for(int i = 0; i < node_count; i++) {
		nodes[i]->project_point_to_line(line_projection_matrix);
	}
}


Node * Element::get_center()
{
	return this->center;
}

double Element::get_distance()
{
	return this->distance;
}

void Element::set_distance(double distance)
{
	this->distance = distance;
}

Element::~Element() {
	delete[] nodes;
	if(center) {
		delete center;
	}
	delete[] closest_elements;
}

int Element::get_element_type(DenseMatrix<int> *els)
{
    if((*els)[ 8*els->get_columns()+0]==0) return M_ELEMENT_LINE2;
    if((*els)[ 9*els->get_columns()+0]==0) return M_ELEMENT_LINE3;
    if((*els)[10*els->get_columns()+0]==0 && (*els)[ 8*els->get_columns()+0]==(*els)[ 9*els->get_columns()+0]) return M_ELEMENT_TRIA3;
    if((*els)[10*els->get_columns()+0]==0) return M_ELEMENT_QUAD4;
    if((*els)[13*els->get_columns()+0] >0 && (*els)[ 8*els->get_columns()+0]==(*els)[ 9*els->get_columns()+0]) return M_ELEMENT_TRIA6;
    if((*els)[13*els->get_columns()+0] >0) return M_ELEMENT_QUAD8;

    return -1;
}

BoundingVolume * Element::get_bound_volume(Node **nodes, int node_count, int bounds_count)
{
    Interval *b = new Interval[bounds_count];

    double min_, max_, value;
    for(int i = 0; i < bounds_count; i++) {
        min_ = Element::get_value_of_fn[i](nodes[0]);
        max_ = min_;
        for(int j = 1; j < node_count; j++) {
            value = Element::get_value_of_fn[i](nodes[j]);
            if(value < min_) {
                min_ = value;
            }
            if(value > max_) {
                max_ = value;
            }
        }
        b[i] = Interval(min_, max_);
    }

    BoundingVolume *bounds = new BoundingVolume(b, bounds_count);
    return bounds;
}

Element * Element::find_neighboor_element(std::map<int, std::vector<Element*> > &adjacentMap, int n1, int n2, int denied_id)
{
    for(size_t i = 0; i < adjacentMap[n1].size(); i++) {
        for(size_t j = 0; j < adjacentMap[n2].size(); j++) {
            int id1 = adjacentMap[n1][i]->get_id();
            int id2 = adjacentMap[n2][j]->get_id();
            if((id1 != denied_id) && (id1 == id2)) {
                return adjacentMap[n1][i];
            }
        }
    }
    return NULL;
}

void Element::add_incident_elements_in_plane(
        std::map<int, std::vector<Element*> > &adjacent_map,
        Element *master,
        Element *slave,
        std::map<int, Element*> *added,
        double *projection_matrix)
{
    bool add = false;
    int mc = master->get_node_count();

    for(int i = 0; i < mc; i++) {
        MCVec3 p = master->get_node(i)->get_plane_projection();
        if(slave->is_point_inside(p)) {
            add = true;
            Node *n = master->get_node(i);
            for(size_t j = 0; j < adjacent_map[n->get_id()].size(); j++) {
                Element *e = adjacent_map[n->get_id()][j];
                if(!added->count(e->get_id())) {
                    e->project_element_to_plane(projection_matrix);
                    added->insert(std::pair<int, Element*>(e->get_id(), e));
                    add_incident_elements_in_plane(adjacent_map, e, slave, added, projection_matrix);
                }
            }
        }
    }

 // if some point is inner we can finish
    if(add) {
        return;
    }
    int sc = slave->get_node_count();
    MCVec3 last;
    for (int i = 0; i < mc; i++) {
        MCVec3 p1 = master->get_node(i)->get_plane_projection();
        MCVec3 p2 = master->get_node((i + 1) % mc)->get_plane_projection();
        MCVec3 u = MCVec3(p2 - p1);
        for (int j = 0; j < sc; j++) {
            MCVec3 q2 = slave->get_node(j)->get_plane_projection();
            MCVec3 v = (q2 - p1);
            if (cross_prod(u, v).z <= 0) {
                MCVec3 prev;
                if (j == 0) {
                    prev = slave->get_node(sc - 1)->get_plane_projection();
                } else {
                    prev = slave->get_node(j - 1)->get_plane_projection();
                }
                MCVec3 next = slave->get_node((j + 1) % sc)->get_plane_projection();
                MCVec3 *i1 = segment_intersect(p1, p2, prev, q2);
                MCVec3 *i2 = segment_intersect(p1, p2, q2, next);
                if (i2 == NULL && i1 == NULL) {
                    break;
                } else {
                    delete i1;
                    delete i2;
                }
                Element *e = find_neighboor_element(
                        adjacent_map,
                        master->get_node(i)->get_id(),
                        master->get_node((i + 1) % mc)->get_id(),
                        master->get_id());
                if (e == NULL) {
                    break;
                }
                if (!added->count(e->get_id())) {
                    e->project_element_to_plane(projection_matrix);
                    added->insert(std::pair<int, Element*>(e->get_id(), e));
                    add_incident_elements_in_plane(adjacent_map, e, slave, added, projection_matrix);
                }
                break;
            }
        }
    }
}

void Element::add_incident_elements_in_line(
        std::map<int, std::vector<Element*> > &adjacent_map,
        Element *master,
        Element *slave,
        std::map<int, Element*> *added,
        double *line_projection_matrix)
{
    bool add = false;
    int mc = master->get_node_count();
    double xm0 = master->get_node(0)->get_plane_projection().x;
    double xmend = master->get_node(master->get_node_count()-1)->get_plane_projection().x;
    double xs0 = slave->get_node(0)->get_plane_projection().x;
    double xsend = slave->get_node(slave->get_node_count()-1)->get_plane_projection().x;
    for(int i = 0; i < mc; i++) {
        double p = master->get_node(i)->get_line_projection().x;
        if(slave->is_point_inside(MCVec3(p, 0.0, 0.0))) {
            add = true;
            Node *n = master->get_node(i);
            for(size_t j = 0; j < adjacent_map[n->get_id()].size(); j++) {
                Element *e = adjacent_map[n->get_id()][j];
                if(!added->count(e->get_id())) {
                    e->project_element_to_line(line_projection_matrix);
                    added->insert(std::pair<int, Element*>(e->get_id(), e));
                    add_incident_elements_in_line(adjacent_map, e, slave, added, line_projection_matrix);
                }
            }
        }
    }
    // if some point is inner we can finish
    if(add) {
        return;
    }
    if ( (std::min(xm0,xmend) <= std::min(xs0,xsend)) && (std::max(xs0,xsend) <= std::max(xm0,xmend)) ) {
    	if (!added->count(master->get_id())) {
    		//master->project_element_to_line(line_projection_matrix);
    		added->insert(std::pair<int, Element*>(master->get_id(), master));
    		//add_incident_elements_in_line(adjacentMap, e, slave, added, line_projection_matrix);
    	}
    }
}

PointList * Element::clip_element_polygons(Element *master, Element *slave)
{
	PointList *master_pl  = new PointList();
	PointList *slave_pl  = new PointList();
	for(int i = 0; i < slave->get_node_count(); i++) {
		slave_pl->push_back(slave->get_node(i)->get_plane_projection());
	}
	for(int i = 0; i < master->get_node_count(); i++) {
		master_pl->push_back(master->get_node(i)->get_plane_projection());
	}
	PointList *clipped;
	switch (slave->get_type()) {
	case M_ELEMENT_LINE2:
	case M_ELEMENT_LINE3: clipped = clip_lines(master_pl, slave_pl); break;
	case M_ELEMENT_TRIA3:
	case M_ELEMENT_TRIA6:
	case M_ELEMENT_QUAD4:
	case M_ELEMENT_QUAD8:
	case M_ELEMENT_QUAD9: clipped = clip_polygons(master_pl, slave_pl); break;
	case M_ELEMENT_UNKNOWN: fprintf(stderr, "Element::clip_element_polygons ... unknown slave element type\n"); exit(1); break;

	}

	return clipped;
}

Element * Element::create_element(int id, Node **nodes, int element_type)
{
	if(element_type == M_ELEMENT_LINE2) {
		return new Element_line2(id, nodes);
	}
	if(element_type == M_ELEMENT_LINE3) {
		return new Element_line3(id, nodes);
	}
	if(element_type == M_ELEMENT_TRIA3) {
		return new Element_tria3(id, nodes);
	}
	if(element_type == M_ELEMENT_TRIA6) {
		return new Element_tria6(id, nodes);
	}
	if(element_type == M_ELEMENT_QUAD4) {
		return new Element_quad4(id, nodes);
	}
	if(element_type == M_ELEMENT_QUAD8) {
		return new Element_quad8(id, nodes);
	}

	return NULL;
}

int Element::get_element_nodes_count(int element_type)
{
	if(element_type == M_ELEMENT_LINE2) {
		return M_LINE2_NODES_COUNT;
	}
	if(element_type == M_ELEMENT_LINE3) {
		return M_LINE3_NODES_COUNT;
	}
	if(element_type == M_ELEMENT_TRIA3) {
		return M_TRIA3_NODES_COUNT;
	}
	if(element_type == M_ELEMENT_TRIA6) {
		return M_TRIA6_NODES_COUNT;
	}
	if(element_type == M_ELEMENT_QUAD4) {
		return M_QUAD4_NODES_COUNT;
	}
	if(element_type == M_ELEMENT_QUAD8) {
		return M_QUAD8_NODES_COUNT;
	}

	return 0;
}

int Element::get_bound_count(int element_type)
{
	if(element_type == M_ELEMENT_LINE2) {
		return M_BOUND_COUNT_2D;
	}
	if(element_type == M_ELEMENT_LINE3) {
		return M_BOUND_COUNT_2D;
	}
	if(element_type == M_ELEMENT_TRIA3) {
		return M_BOUND_COUNT_3D;
	}
	if(element_type == M_ELEMENT_TRIA6) {
		return M_BOUND_COUNT_3D;
	}
	if(element_type == M_ELEMENT_QUAD4) {
		return M_BOUND_COUNT_3D;
	}
	if(element_type == M_ELEMENT_QUAD8) {
		return M_BOUND_COUNT_3D;
	}

	return 0;
}

Value Element::get_value_of_fn[9] = {
		Element::get_value_of_fn_x,
		Element::get_value_of_fn_y,
		Element::get_value_of_fn_y_minus_x,
		Element::get_value_of_fn_y_plus_x,
		Element::get_value_of_fn_z,
		Element::get_value_of_fn_z_minus_x,
		Element::get_value_of_fn_z_plus_x,
		Element::get_value_of_fn_z_minus_y,
		Element::get_value_of_fn_z_plus_y
};

Compare	Element::compare_by_fn[9] = {
		Element::compare_by_fn_x,
		Element::compare_by_fn_y,
		Element::compare_by_fn_y_minus_x,
		Element::compare_by_fn_y_plus_x,
		Element::compare_by_fn_z,
		Element::compare_by_fn_z_minus_x,
		Element::compare_by_fn_z_plus_x,
		Element::compare_by_fn_z_minus_y,
		Element::compare_by_fn_z_plus_y
};

double Element::get_value_of_fn_x(Node *n)
{
	return n->get_coordinates().x;
}

double Element::get_value_of_fn_y(Node *n)
{
	return n->get_coordinates().y;
}

double Element::get_value_of_fn_y_minus_x(Node *n)
{
	return n->get_coordinates().y - n->get_coordinates().x;
}

double Element::get_value_of_fn_y_plus_x(Node *n)
{
	return n->get_coordinates().y + n->get_coordinates().x;
}

double Element::get_value_of_fn_z(Node *n)
{
	return n->get_coordinates().z;
}

double Element::get_value_of_fn_z_minus_x(Node *n)
{
	return n->get_coordinates().z - n->get_coordinates().x;
}

double Element::get_value_of_fn_z_plus_x(Node *n)
{
	return n->get_coordinates().z + n->get_coordinates().x;
}

double Element::get_value_of_fn_z_minus_y(Node *n)
{
	return n->get_coordinates().z - n->get_coordinates().y;
}

double Element::get_value_of_fn_z_plus_y(Node *n)
{
	return n->get_coordinates().z + n->get_coordinates().y;
}


/*
 * Compare two elements by function x
 */
int Element::compare_by_fn_x(const void * element1, const void * element2)
{
	Element **e1 = (Element**)element1;
	Element **e2 = (Element**)element2;
	double d =
			get_value_of_fn[fn_x]((*e1)->get_center()) -
			get_value_of_fn[fn_x]((*e2)->get_center());
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

/*
 * Compare two elements by function y
 */
int Element::compare_by_fn_y(const void * element1, const void * element2)
{
	Element **e1 = (Element**)element1;
	Element **e2 = (Element**)element2;
	double d =
			get_value_of_fn[fn_y]((*e1)->get_center()) -
			get_value_of_fn[fn_y]((*e2)->get_center());
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

/*
 * Compare two elements by function x = y
 */
int Element::compare_by_fn_y_minus_x(const void * element1, const void * element2)
{
	Element **e1 = (Element**)element1;
	Element **e2 = (Element**)element2;
	double d =
			get_value_of_fn[fn_y_minus_x]((*e1)->get_center()) -
			get_value_of_fn[fn_y_minus_x]((*e2)->get_center());
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

/*
 * Compare two elements by function -x = y
 */
int Element::compare_by_fn_y_plus_x(const void * element1, const void * element2)
{
	Element **e1 = (Element**)element1;
	Element **e2 = (Element**)element2;
	double d =
			get_value_of_fn[fn_y_plus_x]((*e1)->get_center()) -
			get_value_of_fn[fn_y_plus_x]((*e2)->get_center());
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

/*
 * Compare two elements by function z
 */
int Element::compare_by_fn_z(const void * element1, const void * element2)
{
	Element **e1 = (Element**)element1;
	Element **e2 = (Element**)element2;
	double d =
			get_value_of_fn[fn_z]((*e1)->get_center()) -
			get_value_of_fn[fn_z]((*e2)->get_center());
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

/*
 * Compare two elements by function x = z
 */
int Element::compare_by_fn_z_minus_x(const void * element1, const void * element2)
{
	Element **e1 = (Element**)element1;
	Element **e2 = (Element**)element2;
	double d =
			get_value_of_fn[fn_z_minus_x]((*e1)->get_center()) -
			get_value_of_fn[fn_z_minus_x]((*e2)->get_center());
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

/*
 * Compare two elements by function -x = z
 */
int Element::compare_by_fn_z_plus_x(const void * element1, const void * element2)
{
	Element **e1 = (Element**)element1;
	Element **e2 = (Element**)element2;
	double d =
			get_value_of_fn[fn_z_plus_x]((*e1)->get_center()) -
			get_value_of_fn[fn_z_plus_x]((*e2)->get_center());
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

/*
 * Compare two elements by function z = y
 */
int Element::compare_by_fn_z_minus_y(const void * element1, const void * element2)
{
	Element **e1 = (Element**)element1;
	Element **e2 = (Element**)element2;
	double d =
			get_value_of_fn[fn_z_minux_y]((*e1)->get_center()) -
			get_value_of_fn[fn_z_minux_y]((*e2)->get_center());
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}

/*
 * Compare two elements by function z = -y
 */
int Element::compare_by_fn_z_plus_y(const void * element1, const void * element2)
{
	Element **e1 = (Element**)element1;
	Element **e2 = (Element**)element2;
	double d =
			get_value_of_fn[fn_z_plus_y]((*e1)->get_center()) -
			get_value_of_fn[fn_z_plus_y]((*e2)->get_center());
	if(d <= 0) {
		return -1;
	} else {
		return 1;
	}
}
