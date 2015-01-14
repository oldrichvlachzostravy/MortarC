/*
 * BoundingVolume.h
 *
 *  Created on: Oct 29, 2014
 *      Author: bullseye
 */

#ifndef BOUNDINGVOLUME_H_
#define BOUNDINGVOLUME_H_

#include <assert.h>

#define EPSILON 0.01

/**
 * Structure that represents the interval on the real axis (\f$ I=\langle a,b \rangle=\{x\in\mathbb{R}|\, a \leq x \leq b\}\f$)
 */
struct Interval
{
        /**
         * default constructor \f$ a=b=0 \f$
         */
        Interval(): start(0), end(0) { };
        /**
         * constructor
         * @param s = center
         * @param e = (a+b)/2
         */
        Interval(double s, double e): start(s - EPSILON), end(e + EPSILON) {
            assert(start < end);
        };

        double start;
        double end;
        /**
         * @return true if this is the subset of interval, \f$ this \subseteq interval \f$, false otherwise
         * @param interval given interval
         */
        bool is_overlapped(Interval interval) {
            if(interval.start > end) return false;
            if(interval.end < start) return false;
            return true;
        }

        /**
         * @return the size of this
         */
        double get_interval_size() { return end - start; }
};


/**
 * The 3D \f$ k--DOP\f$
 */
class BoundingVolume
{
    public:
        BoundingVolume(Interval*, int);
        ~BoundingVolume();

        Interval * get_bounds();
        Interval get_bound(int);
        int get_bounds_count();
        double get_biggest_interval_size();
        int get_biggest_interval_index();

        bool is_overlapped(BoundingVolume*);

    private:
        int bounds_count;
        Interval *bounds;
};


#endif /* BOUNDINGVOLUME_H_ */
