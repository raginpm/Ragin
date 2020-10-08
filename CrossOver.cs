using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ConsoleApp42
{
    class CrossOver
    {
        public int[] Cross(float[] MA1, float[] MA2)
        {
            int[] crs = new int[MA1.Length];
            Array.Clear(crs, 0, crs.Length);
            for (int i = 1; i < MA1.Length; i++)
            {
                crs[i] = 0;
                if (MA1[i - 1] != 0 && MA2[i - 1] != 0)
                {
                    if ((MA1[i - 1] < MA2[i - 1] && MA1[i] > MA2[i]) || (MA1[i - 1] == MA2[i - 1] && MA1[i] > MA2[i]))
                    {
                        crs[i] = 1;
                    }
                }
            }
            return crs;
        }
    }
}
