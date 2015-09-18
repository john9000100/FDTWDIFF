package com.util;

public class DotProductDistance implements DistanceFunction
{

   public DotProductDistance()
   {
      
   }
   
   
   public double calcDistance(double[] vector1, double[] vector2)
   {
      if (vector1.length != vector2.length)
         throw new InternalError("ERROR:  cannot calculate the distance "
                                    + "between vectors of different sizes.");

      double dotProduct = 0.0;
      for (int x=0; x<vector1.length; x++)
         dotProduct += vector1[x]*vector2[x];

      //negate and shift the dot product
      dotProduct = -(dotProduct - 1);

      return dotProduct;
   }

}
