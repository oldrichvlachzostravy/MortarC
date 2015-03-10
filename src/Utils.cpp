#include "Utils.h"


MCVec3 transform(MCVec3 p, double *matrix)
{
	MCVec3 result;
	result.x = p.x * matrix[0] + p.y * matrix[3] + p.z * matrix[6];
	result.y = p.x * matrix[1] + p.y * matrix[4] + p.z * matrix[7];
	result.z = p.x * matrix[2] + p.y * matrix[5] + p.z * matrix[8];
	return result;
}

MCVec2 transform(MCVec2 p, double *matrix)
{
	MCVec2 result;
	result.x = p.x * matrix[0] + p.y * matrix[2];
	result.y = p.x * matrix[1] + p.y * matrix[3];
	return result;
}

MCVec3 back_transform(MCVec3 p, double *matrix)
{
	MCVec3 result;
	result.x = p.x * matrix[0] + p.y * matrix[1] + p.z * matrix[2];
	result.y = p.x * matrix[3] + p.y * matrix[4] + p.z * matrix[5];
	result.z = p.x * matrix[6] + p.y * matrix[7] + p.z * matrix[8];
	return result;
}

MCVec2 back_transform(MCVec2 p, double *matrix)
{
	MCVec2 result;
	result.x = p.x * matrix[0] + p.y * matrix[1];
	result.y = p.x * matrix[2] + p.y * matrix[3];
	return result;
}
MCVec3 * segment_intersect(MCVec3 p0, MCVec3 p1, MCVec3 q0, MCVec3 q1)
{
	MCVec3 u = p1 - p0;
	MCVec3 v = q1 - q0;
	MCVec3 w = p0 - q0;
	const double e = 0.0001;

	double scale = u.x * v.y - u.y * v.x;
	if(fabs(scale) < e * u.length() * v.length()) { // lines are parallel
		return NULL;
	}

	double s = (w.y * v.x - w.x * v.y) / scale;
	if(s > 1 + e || s < 0 - e) {
		return NULL;
	}

	double t = (w.y * u.x - w.x * u.y) / scale;
	if(t > 1 + e || t < 0 - e) {
		return NULL;
	}

	MCVec3 *intersect = new MCVec3(q0 + v * t);
	return intersect;
}

MCVec3 * line_intersect(MCVec3 p0, MCVec3 p1, MCVec3 q0, MCVec3 q1)
{
	MCVec3 u = p1 - p0;
	MCVec3 v = q1 - q0;
	MCVec3 w = p0 - q0;
	const double e = 0.0001;

	double scale = u.x * v.y - u.y * v.x;
	if(fabs(scale) < e * u.length() * v.length()) { // lines are parallel
		return NULL;
	}

	double t = (w.y * u.x - w.x * u.y) / scale;

	MCVec3 *intersect = new MCVec3(q0 + v * t);
	return intersect;
}

MCVec3 * line_plane_intersect(MCVec3 p0, MCVec3 p1, MCVec3 v0, MCVec3 v1, MCVec3 v2)
{
	double e = 0.001;
	MCVec3 u = v1 - v0;
	MCVec3 v = v2 - v0;
	MCVec3 p = p1 - p0;
	MCVec3 n = cross_prod(u, v);

	double s = dot_prod(n, p);
	if(fabs(s) < e) {
		//TODO: test if line is in plane
		return NULL;
	}

	double t = dot_prod(n, (v0 - p0)) / s;
	MCVec3 *intersect = new MCVec3(p0 + p * t);
	return intersect;
}

MCVec3 get_barycentric(MCVec3 u, MCVec3 v, MCVec3 w)
{
	// from http://geomalgorithms.com/a04-_planes.html#Barycentric-Coordinates
	double uv = dot_prod(u, v);
	double wu = dot_prod(w, u);
	double wv = dot_prod(w, v);
	double uu = dot_prod(u, u);
	double vv = dot_prod(v, v);
	double scale = uv * uv - uu * vv;
	double s = (uv * wv - vv * wu) / scale;
	double t = (uv * wu - uu * wv) / scale;
	return MCVec3(1.0 - s - t, s, t);
}

MCVec3 get_reference(MCVec3 u, MCVec3 v, MCVec3 w)
{
	double uv = dot_prod(u, v);
	double wu = dot_prod(w, u);
	double wv = dot_prod(w, v);
	double uu = dot_prod(u, u);
	double vv = dot_prod(v, v);
	double scale = uv * uv - uu * vv;
	double s = (uv * wv - vv * wu) / scale;
	double t = (uv * wu - uu * wv) / scale;
	return MCVec3(s, t, 0.0);
}

PointList* clip_lines(PointList *master, PointList *slave){
	MCVec3* left, *right;
	if (slave->begin()->x < slave->rbegin()->x) {
		if (master->begin()->x < master->rbegin()->x) {
			left  = (slave-> begin()->x < master-> begin()->x)? &(master-> begin().operator *()) :  &( slave-> begin().operator *());
			right = (slave->rbegin()->x < master->rbegin()->x)? &( slave->rbegin().operator *()) :  &(master->rbegin().operator *());
		} else {
			left  = (slave-> begin()->x < master->rbegin()->x)? &(master->rbegin().operator *()) :  &( slave-> begin().operator *());
			right = (slave->rbegin()->x < master-> begin()->x)? &( slave->rbegin().operator *()) :  &(master-> begin().operator *());
		}
	} else {
		if (master->begin()->x < master->rbegin()->x) {
			left  = (slave->rbegin()->x < master-> begin()->x)? &(master-> begin().operator *()) :  &( slave->rbegin().operator *());
			right = (slave-> begin()->x < master->rbegin()->x)? &( slave-> begin().operator *()) :  &(master->rbegin().operator *());
		} else {
			left  = (slave->rbegin()->x < master->rbegin()->x)? &(master->rbegin().operator *()) :  &( slave-> begin().operator *());
			right = (slave-> begin()->x < master-> begin()->x)? &( slave-> begin().operator *()) :  &(master-> begin().operator *());
		}
	}
	PointList *clipped = new PointList();
	clipped->push_back(*left);
	clipped->push_back(*right);
//	printf("slave  [%f, %f]--[%f, %f]\n",
//			slave->begin().operator *().x, slave->begin().operator *().y,
//			slave->begin().operator ++().operator *().x, slave->begin().operator ++().operator *().y);
//	printf("master [%f, %f]--[%f, %f]\n",
//				master->begin().operator *().x, master->begin().operator *().y,
//				master->begin().operator ++().operator *().x, master->begin().operator ++().operator *().y);
	return clipped;
}

PointList* clip_polygons(PointList *master, PointList *slave){
	MCVec3 last, start;
	int previous, first;
	//int nodes = master->size();
	PointList *clipped = new PointList(slave->begin(), slave->end());

	PointList::iterator clipped_it;

	for(PointList::iterator master_it = master->begin(); master_it !=master->end(); master_it++) {
		MCVec3 p2 = *master_it;
		PointList::iterator master_it_next_or_begin = master_it;
		master_it_next_or_begin++;
		if (master_it_next_or_begin == master->end()) {
			master_it_next_or_begin = master->begin();
		}
		MCVec3 p1 = *master_it_next_or_begin;
		MCVec3 u = MCVec3(p2 - p1);

		previous = -1;
		first = -1;
		for(clipped_it = clipped->begin(); clipped_it != clipped->end(); clipped_it++) {
			if(previous == -1) {
				start = (*clipped_it);
			}

			MCVec3 q2 = (*clipped_it);
			MCVec3 v = MCVec3(q2 - p1);
			if(cross_prod(u, v).z >= 0) {
				if(previous == 0) {
					MCVec3 *p = line_intersect(p1, p2, last, q2);
					if (p) {
						clipped->insert(clipped_it, *p);
						delete p;
					}
				}
				last = (*clipped_it);
				previous = 1;
			} else {
				if(previous == 1) {
					MCVec3 *p = line_intersect(p1, p2, last, q2);
					if(p) {
						clipped->insert(clipped_it, *p);
						delete p;
					}
				}
				last = (*clipped_it);
				clipped->erase(clipped_it--);
				previous = 0;
			}
			if(first == -1) {
				first = previous;
			}
		}
		if(previous != first) {
			MCVec3 *p = line_intersect(p1, p2, last, start);
			if (p) {
				clipped->insert(clipped_it, *p);
				delete p;
			}
		}
	}
	// remove points that are too close
	if (clipped->size() > 0) {
		double e = 0.000001; // minimal points distance
		MCVec3 l, f;
		clipped_it = clipped->begin();
		l = f = *clipped_it;
		for (clipped_it++; clipped_it != clipped->end(); clipped_it++) {
			if (((l.x < clipped_it->x + e) && (l.x > clipped_it->x - e)) && ((l.y < clipped_it->y + e) && (l.y > clipped_it->y - e))) {
				clipped->erase(clipped_it--);
			}
			l = *clipped_it;
		}
		if (clipped->size() > 0 && ((l.x < f.x + e) && (l.x > f.x - e)) && ((l.y < f.y + e) && (l.y > f.y - e))) {
			clipped->erase(clipped->begin());
		}
	}
	return clipped;
}

bool legalize_edge(MCVec3 pi, MCVec3 pj, MCVec3 pk, MCVec3 pl)
{
	MCVec3 vecij = MCVec3((pj - pi));
	MCVec3 midij = MCVec3(pi + vecij / 2);
	MCVec3 tecij = MCVec3(vecij.y, -vecij.x, vecij.z);

	MCVec3 vecik = MCVec3((pk - pi));
	MCVec3 midik = MCVec3(pi + vecik / 2);
	MCVec3 tecik = MCVec3(vecik.y, -vecik.x, vecik.z);

	MCVec3 *center = line_intersect(midij, midij + tecij, midik, midik + tecik);

	double rr = ((*center) - pi).length();
	if(rr <= ((*center) - pl).length()) {
		return true;
	} else {
		return false;
	}
}

void legalize_edge(TDescription *t1, TDescription *t2, int side)
{
	MCVec3 vecij = MCVec3((t1->p[1] - t1->p[0]));
	MCVec3 midij = MCVec3(t1->p[0] + vecij / 2);
	MCVec3 tecij = MCVec3(vecij.y, -vecij.x, vecij.z);

	MCVec3 vecik = MCVec3((t1->p[2] - t1->p[0]));
	MCVec3 midik = MCVec3(t1->p[0] + vecik / 2);
	MCVec3 tecik = MCVec3(vecik.y, -vecik.x, vecik.z);

	MCVec3 *center = line_intersect(midij, midij + tecij, midik, midik + tecik);

	double rr;
	if(center) {
		rr = ((*center) - t1->p[0]).length();
	}
	if(center && rr > ((*center) - t2->p[t1->t[side]->first]).length()) {
		int idx2 = t1->t[side]->first;
		int idx1 = t2->t[idx2]->first;
		t2->p[(idx2 + 2) % 3] = t1->p[idx1];
		t1->p[(idx1 + 2) % 3] = t2->p[idx2];
		TDescription *prev1 = NULL;
		TDescription *prev2 = NULL;
		if(t1->t[(idx1 + 1) % 3] != NULL) {
			prev1 = t1->t[(idx1 + 1) % 3]->second;
			prev1->t[t1->t[(idx1 + 1) % 3]->first]->second = t2;
		}
		if (t1->t[(idx1 + 2) % 3] != NULL) {
			prev2 = t1->t[(idx1 + 2) % 3]->second;
			prev2->t[t1->t[(idx1 + 2) % 3]->first]->second = t1;
		}
		delete t2->t[idx2];
		delete t1->t[idx1];
		t2->t[idx2] = t1->t[(idx1 + 1) % 3];
		t1->t[idx1] = t2->t[(idx2 + 1) % 3];
		t2->t[(idx2 + 1) % 3] = new std::pair<int, TDescription*>((idx1 + 1) % 3, t1);
		t1->t[(idx1 + 1) % 3] = new std::pair<int, TDescription*>((idx2 + 1) % 3, t2);

		if(prev1 != NULL) {
			legalize_edge(prev1, t2, t2->t[idx2]->first);
		}
		if(prev2 != NULL) {
			legalize_edge(prev2, t1, t1->t[(idx1 + 2) % 3]->first);
		}
	}
	if(center) {
		delete center;
	}
}

std::vector<TDescription*> triangulate(PointList* polygon)
{
	//mexPrintf("triangulation.....(");
	std::vector<TDescription*> triangles = std::vector<TDescription*>();

	PointList::iterator it;
	/*it = polygon->begin();
	while (it != polygon->end()) {
		mexPrintf("[%f, %f] ", it->x, it->y);
		it++;
	}
	mexPrintf("---> ");*/

	it = polygon->begin();
	TDescription *t = new TDescription();
	t->p[0] = *it++;
	t->p[1] = *it++;
	t->p[2] = *it++;
	triangles.push_back(t);

	int i = 0;
	while (it != polygon->end()) {
		//mexPrintf("[%f, %f] ", it->x, it->y);
		TDescription *t = new TDescription();
		t->p[0] = triangles[i]->p[0];
		t->p[1] = triangles[i]->p[2];
		t->p[2] = *it++;
		triangles[i]->t[1] = new std::pair<int, TDescription*>(2, t);
		t->t[2] = new std::pair<int, TDescription*>(1, triangles[i]);
		triangles.push_back(t);
		i++;
		legalize_edge(triangles[i - 1], triangles[i], 1);
	}

	return triangles;

}

void dense_matrix_inverse(double * a, int n)
{
	LAPACKINT nn = n;
	LAPACKINT *ipiv = new LAPACKINT[nn+1];
	LAPACKINT lwork = nn*nn;
	double *work = new double[lwork];
	LAPACKINT info;

	dgetrf_(&nn,&nn,a,&nn,ipiv,&info);
	dgetri_(&nn,a,&nn,ipiv,work,&lwork,&info);

	delete ipiv;
	delete work;
}

void dense_matrix_multiply(double * a, double * b, double * c, int m, int n, int o)
{
	const char trans = 'N';
	double alpha = 1.0;
	double beta  = 0.0;

	if (o==0) o = m;
	if (n==0) n = o;

	dgemm_(&trans, &trans, &m, &o, &n, &alpha, a, &m, b, &n, &beta, c, &m);
}

void dense_matrix_solve(double * a, double * b, int n, int o)
{
	LAPACKINT nn = n;
	LAPACKINT oo = o;
	LAPACKINT *ipiv = new LAPACKINT[nn+1];
	LAPACKINT info;

	dgesv_(&nn, &oo, a, &nn, ipiv, b, &nn, &info);
	delete[] ipiv;
}

void dense_matrix_transpose(double * a, int m, int n)
{
	for (int i=0; i < m; i++)
	{
		for (int j=i+1; j < n; j++)
		{
			double tmp = a[j*m+i];
			a[j*m+i] = a[i*n+j];
			a[i*n+j] = tmp;
		}
	}
}

void mex_print_full_matrix(double * a, int m, int n)
{
//	for (int i = 0; i < m; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			mexPrintf(" %.3f", a[j*m+i]);
//		}
//		mexPrintf("\n");
//	}
}

void print_sparse_matrix(std::map<int,std::map<int,double> > &m, const char *name)
{
	std::cout << name << std::endl;
	for (std::map<int,std::map<int,double> >::iterator rit = m.begin(); rit != m.end(); ++rit) {
		for (std::map<int,double>::iterator cit = rit->second.begin(); cit != rit->second.end(); ++cit) {
			std::cout << std::setw(4) << rit->first << ", " << std::setw(4) << cit->first << ", " << std::setw(7) << cit->second << std::endl;
		}
	}
}
