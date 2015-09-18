/*
 * FastLightDTW.java   Jul 14, 2004
 *
 * Copyright (c) 2004 Stan Salvador
 * stansalvador@hotmail.com
 */

package com.dtw;

import com.timeseries.TimeSeries;
import com.timeseries.PAA;
import com.util.DistanceFunction;


public class FastLightDTW
{
   // CONSTANTS
   final static int DEFAULT_SEARCH_RADIUS = 1;


   public static double getWarpDistBetween(TimeSeries tsI, TimeSeries tsJ, DistanceFunction distFn, boolean flexible)
   {
      return fastLightDTW(tsI, tsJ, DEFAULT_SEARCH_RADIUS, distFn, flexible).getDistance();
   }


   public static double getWarpDistBetween(TimeSeries tsI, TimeSeries tsJ, int searchRadius, DistanceFunction distFn, boolean flexible)
   {
      return fastLightDTW(tsI, tsJ, searchRadius, distFn, flexible).getDistance();
   }


   public static WarpPath getWarpPathBetween(TimeSeries tsI, TimeSeries tsJ, DistanceFunction distFn, boolean flexible)
   {
      return fastLightDTW(tsI, tsJ, DEFAULT_SEARCH_RADIUS, distFn, flexible).getPath();
   }


   public static WarpPath getWarpPathBetween(TimeSeries tsI, TimeSeries tsJ, int searchRadius, DistanceFunction distFn, boolean flexible)
   {
      return fastLightDTW(tsI, tsJ, searchRadius, distFn, flexible).getPath();
   }


   public static TimeWarpInfo getWarpInfoBetween(TimeSeries tsI, TimeSeries tsJ, int searchRadius, DistanceFunction distFn, boolean flexible)
   {
      return fastLightDTW(tsI, tsJ, searchRadius, distFn, flexible);
   }


   private static TimeWarpInfo fastLightDTW(TimeSeries tsI, TimeSeries tsJ, int searchRadius, DistanceFunction distFn, boolean flexible)
   {
      if (searchRadius < 0)
         searchRadius = 0;

      final int minTSsize = searchRadius+2;

      if ( (tsI.size() <= minTSsize) || (tsJ.size()<=minTSsize) )
      {
         // Perform full Dynamic Time Warping.
         return LightDTW.getWarpInfoBetween(tsI, tsJ, distFn, flexible);
      }
      else
      {
         final double resolutionFactor = 2.0;

         final PAA shrunkI = new PAA(tsI, (int)(tsI.size()/resolutionFactor));
         final PAA shrunkJ = new PAA(tsJ, (int)(tsJ.size()/resolutionFactor));

 
          // Determine the search window that constrains the area of the cost matrix that will be evaluated based on
          //    the warp path found at the previous resolution (smaller time series).
	WarpPath shrunkWarpPath = FastLightDTW.getWarpPathBetween(shrunkI, shrunkJ, searchRadius, distFn, flexible);

        final SearchWindow window = new ExpandedResWindow(tsI, tsJ, shrunkI, shrunkJ,
                                                            shrunkWarpPath,
                                                            searchRadius);



         // Find the optimal warp path through this search window constraint.
         return LightDTW.getWarpInfoBetween(tsI, tsJ, window, distFn, flexible);
      }  // end if
   }  // end recFastLightDTW(...)

}  // end class FastLightDTW
