/*
 * FastDtwTest.java   Jul 14, 2004
 *
 * Copyright (c) 2004 Stan Salvador
 * stansalvador@hotmail.com
 */

package com;

import com.timeseries.TimeSeries;
import com.timeseries.TimeSeriesPoint;
import com.util.DistanceFunction;
import com.util.DistanceFunctionFactory;
import com.dtw.TimeWarpInfo;
import com.dtw.MATLABTimeWarpInfo;
import java.util.ArrayList;


/**
 * This class contains a main method that executes the FastDTW algorithm on two
 * time series with a specified radius.
 *
 * @author Stan Salvador, stansalvador@hotmail.com
 * @since Jul 14, 2004
 */
public class FastDtwTest
{
   /**
    * This main method executes the FastDTW algorithm on two time series with a
    * specified radius. The time series arguments are file names for files that
    * contain one measurement per line (time measurements are an optional value
    * in the first column). After calculating the warp path, the warp
    * path distance will be printed to standard output, followed by the path
    * in the format "(0,0),(1,0),(2,1)..." where each pair of numbers in
    * parenthesis are indexes of the first and second time series that are
    * linked in the warp path
    *
    * Methods for calling DTW and Fast DTW from MATLAB are below the main method.
    * @param args  command line arguments (see method comments)
    */
      public static void main(String[] args)
      {
         if (args.length!=3 && args.length!=4)
         {
            System.out.println("USAGE:  java FastDtwTest timeSeries1 timeSeries2 radius [EuclideanDistance|ManhattanDistance|BinaryDistance]");
            System.exit(1);
         }
         else
         {
            final TimeSeries tsI = new TimeSeries(args[0], false, false, ',');
            final TimeSeries tsJ = new TimeSeries(args[1], false, false, ',');
            
            final DistanceFunction distFn;
            if (args.length < 4)
            {
               distFn = DistanceFunctionFactory.getDistFnByName("EuclideanDistance"); 
            }
            else
            {
               distFn = DistanceFunctionFactory.getDistFnByName(args[3]);
            }   // end if
            
            final TimeWarpInfo info = com.dtw.FastDTW.getWarpInfoBetween(tsI, tsJ, Integer.parseInt(args[2]), distFn);

            System.out.println("Warp Distance: " + info.getDistance());
            System.out.println("Warp Path:     " + info.getPath());
         }  // end if

      }  // end main()


/* The following methods are for calling from MATLAB */

/** Run fast DTW.
* @param mlq        (# of frames) x (# of dimensions) query sequence
* @param mld        (# of frames) x (# of dimensions) database sequence  
* @param distFunc   Distance Function. See src/com/util for the DistanceFunction classes.
*/
      public static MATLABTimeWarpInfo fastDTW(double[][] mlq, double[][] mld, int radius, String distFunc) {
	 return fastDTW(mlq, mld, radius, distFunc, 0, 0, 0, "Dot", false);
      }

/** Run fast DTW.
* @param mlq            (# of frames) x (# of dimensions) query sequence
* @param mld            (# of frames) x (# of dimensions) database sequence  
* @param radius         window radius for Fast DTW.
* @param distFunc       Distance Function. See src/com/util for the DistanceFunction classes.
* @param vcost          a constant cost associated with moving from (i, j) to (i-1, j), duplicating a database frame.
* @param hcost          a constant cost associated with moving from (i, j) to (i, j-1), duplicating a query frame. 
* @param motionWeight   a multiplicative factor for the motion weights, which are penalties calculated during warping, based on query and database dynamics.
* @param motionDist     either "Dot" or "Euc", indicating the distance measure to use for calculating the motion weights.
* @param normalize      whether to normalize the warp path score by the warp path length. This is useful in flexible matching, where all
*                       the query frames must be used, but not all the database frames, to avoid the inherent bias against horizontal movements, which would increase the warp path length.
*                       The method does not compute the optimal normalized path.
*/
      public static MATLABTimeWarpInfo fastDTW(double[][] mlq, double[][] mld, int radius, String distFunc, double vcost, double hcost, double motionWeight, String motionDist, boolean normalize) 
      {

	 ArrayList<TimeSeriesPoint> q = new ArrayList<TimeSeriesPoint>();
	 ArrayList<TimeSeriesPoint> d = new ArrayList<TimeSeriesPoint>();	

	 for (int i = 0; i < mlq.length; i++) {
		q.add(new TimeSeriesPoint(mlq[i]));
	 }
	 for (int i = 0; i < mld.length; i++) {
		d.add(new TimeSeriesPoint(mld[i]));
	 }

            final TimeSeries tsI = new TimeSeries(q);
            final TimeSeries tsJ = new TimeSeries(d);
            
            final DistanceFunction distFn;
            
            distFn = DistanceFunctionFactory.getDistFnByName(distFunc);

            final TimeWarpInfo info = com.dtw.FastDTW.getWarpInfoBetween(tsI, tsJ, radius, distFn, vcost, hcost, motionWeight, motionDist, normalize);

	return new MATLABTimeWarpInfo(info);
           
      }  

/** Run DTW.
* @param mlq        (# of frames) x (# of dimensions) query sequence
* @param mld        (# of frames) x (# of dimensions) database sequence  
* @param distFunc   Distance Function. See src/com/util for the DistanceFunction classes.
*/
      public static MATLABTimeWarpInfo DTW(double[][] mlq, double[][] mld, String distFunc) {
	 return DTW(mlq, mld, distFunc, 0, 0, 0, "Dot", false);
      }

/** Run DTW.
* @param mlq           (# of frames) x (# of dimensions) query sequence
* @param mld           (# of frames) x (# of dimensions) database sequence  
* @param distFunc       Distance Function. See src/com/util for the DistanceFunction classes.
* @param vcost          a constant cost associated with moving from (i, j) to (i-1, j), duplicating a database frame.
* @param hcost          a constant cost associated with moving from (i, j) to (i, j-1), duplicating a query frame. 
* @param motionWeight   a multiplicative factor for the motion weights, which are penalties calculated during warping, based on query and database dynamics.
* @param motionDist     either "Dot" or "Euc", indicating the distance measure to use for calculating the motion weights.
* @param normalize      whether to normalize the warp path score by the warp path length. This is useful in flexible matching, where all
*                       the query frames must be used, but not all the database frames, to avoid the inherent bias against horizontal movements, which would increase the warp path length.
*                       The method does not compute the optimal normalized path.
*/
      public static MATLABTimeWarpInfo DTW(double[][] mlq, double[][] mld, String distFunc, double vcost, double hcost, double motionWeight, String motionDist, boolean normalize) {
	 ArrayList<TimeSeriesPoint> q = new ArrayList<TimeSeriesPoint>();
	 ArrayList<TimeSeriesPoint> d = new ArrayList<TimeSeriesPoint>();	

	 for (int i = 0; i < mlq.length; i++) {
		q.add(new TimeSeriesPoint(mlq[i]));
	 }
	 for (int i = 0; i < mld.length; i++) {
		d.add(new TimeSeriesPoint(mld[i]));
	 }

        final TimeSeries tsI = new TimeSeries(q);
        final TimeSeries tsJ = new TimeSeries(d);
            
        final DistanceFunction distFn;
            
        distFn = DistanceFunctionFactory.getDistFnByName(distFunc); 

	final TimeWarpInfo info = com.dtw.DTW.getWarpInfoBetween(tsI, tsJ, distFn, vcost, hcost, motionWeight, motionDist, normalize);

	return new MATLABTimeWarpInfo(info);
      }

     /* Method for calling full flexible DTW; this only uses the endpoint modification and is much faster than the version with motion weights. */
     public static MATLABTimeWarpInfo flexibleDTW(double[][] mlq, double[][] mld, String distFunc) {
	 ArrayList<TimeSeriesPoint> q = new ArrayList<TimeSeriesPoint>();
	 ArrayList<TimeSeriesPoint> d = new ArrayList<TimeSeriesPoint>();	

	 for (int i = 0; i < mlq.length; i++) {
		q.add(new TimeSeriesPoint(mlq[i]));
	 }

	 for (int i = 0; i < mld.length; i++) {
		d.add(new TimeSeriesPoint(mld[i]));
	 }

         final TimeSeries tsI = new TimeSeries(q);
         final TimeSeries tsJ = new TimeSeries(d);
            
         final DistanceFunction distFn;
            
         distFn = DistanceFunctionFactory.getDistFnByName(distFunc); 

         final TimeWarpInfo info = com.dtw.LightDTW.getWarpInfoBetween(tsI, tsJ, distFn, true);

	 return new MATLABTimeWarpInfo(info);
      }

      /* Method for calling fast flexible DTW; this only uses the endpoint modification and is much faster than the version with motion weights. */
      public static MATLABTimeWarpInfo fastFlexibleDTW(double[][] mlq, double[][] mld, int radius, String distFunc)
      {
	 ArrayList<TimeSeriesPoint> q = new ArrayList<TimeSeriesPoint>();
	 ArrayList<TimeSeriesPoint> d = new ArrayList<TimeSeriesPoint>();	

	 for (int i = 0; i < mlq.length; i++) {
		q.add(new TimeSeriesPoint(mlq[i]));
	 }

	 for (int i = 0; i < mld.length; i++) {
		d.add(new TimeSeriesPoint(mld[i]));
	 }

         final TimeSeries tsI = new TimeSeries(q);
         final TimeSeries tsJ = new TimeSeries(d);
            
         final DistanceFunction distFn;
            
         distFn = DistanceFunctionFactory.getDistFnByName(distFunc); 

         final TimeWarpInfo info = com.dtw.FastLightDTW.getWarpInfoBetween(tsI, tsJ, radius, distFn, true);

	 return new MATLABTimeWarpInfo(info);

      }  

      /* Method for calling inflexible DTW (the original version, with endpoint constraints) */
      public static MATLABTimeWarpInfo inflexibleDTW(double[][] mlq, double[][] mld, String distFunc) {
	 ArrayList<TimeSeriesPoint> q = new ArrayList<TimeSeriesPoint>();
	 ArrayList<TimeSeriesPoint> d = new ArrayList<TimeSeriesPoint>();	

	 for (int i = 0; i < mlq.length; i++) {
		q.add(new TimeSeriesPoint(mlq[i]));
	 }

	 for (int i = 0; i < mld.length; i++) {
		d.add(new TimeSeriesPoint(mld[i]));
	 }

         final TimeSeries tsI = new TimeSeries(q);
         final TimeSeries tsJ = new TimeSeries(d);
            
         final DistanceFunction distFn;
            
         distFn = DistanceFunctionFactory.getDistFnByName(distFunc); 

         final TimeWarpInfo info = com.dtw.LightDTW.getWarpInfoBetween(tsI, tsJ, distFn, false);

	 return new MATLABTimeWarpInfo(info);
      }

      /* Method for calling fast inflexible DTW (the original version, with endpoint constraints) */
      public static MATLABTimeWarpInfo fastInflexibleDTW(double[][] mlq, double[][] mld, int radius, String distFunc)
      {
	 ArrayList<TimeSeriesPoint> q = new ArrayList<TimeSeriesPoint>();
	 ArrayList<TimeSeriesPoint> d = new ArrayList<TimeSeriesPoint>();	

	 for (int i = 0; i < mlq.length; i++) {
		q.add(new TimeSeriesPoint(mlq[i]));
	 }

	 for (int i = 0; i < mld.length; i++) {
		d.add(new TimeSeriesPoint(mld[i]));
	 }

         final TimeSeries tsI = new TimeSeries(q);
         final TimeSeries tsJ = new TimeSeries(d);
            
         final DistanceFunction distFn;
            
         distFn = DistanceFunctionFactory.getDistFnByName(distFunc); 

         final TimeWarpInfo info = com.dtw.FastLightDTW.getWarpInfoBetween(tsI, tsJ, radius, distFn, false);

	 return new MATLABTimeWarpInfo(info);

      }  

}  // end class FastDtwTest
