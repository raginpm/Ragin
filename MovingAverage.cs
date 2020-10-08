using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ConsoleApp42
{
    class MovingAverage
    {

        public float[] Calc(float[] Close,int Roll)
        {
            float[] MA = new float[Close.Length];
            for (int i = Roll - 1; i < Close.Length; i++)
            {
                float[] sum = new float[1] { 0 };

                for (int j = 0; j < Roll; j++)
                {
                    sum[0] = sum[0] + Close[i - j];
                }
                MA[i] = sum[0] / Roll;
            }
            return MA;
        }
    }
}
