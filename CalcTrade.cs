using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Globalization;
using System.Collections;

namespace ConsoleApp42
{
    class CalcTrade
    {
        List<int> EntrylistCandle = new List<int>();
        List<int> ExitlistCandle = new List<int>();
        List<String> Strategy = new List<String>();
        List<int> Type = new List<int>();
        public void Report(float[] Close,string[] Date)
        {
            foreach (var item in Entry.EntryList.ToArray())
            {
                var i = item.Candle;
                var j = item.SomeValue;
                var k = item.Type;
                EntrylistCandle.Add(i);
                Type.Add(j);
                Strategy.Add(k);
            }
            foreach (var item in Exit.ExitList.ToArray())
            {
                var i = item.Candle;
                ExitlistCandle.Add(i);
            }
            int[] En = EntrylistCandle.ToArray();
            int[] Ex = ExitlistCandle.ToArray();
            var file = @"C:\Users\ragin.pm\Documents\GoldReport1.csv";
            using (var stream = File.CreateText(file))
            {
                for (int i = 0; i < En.Length; i++)
                {
                    if (i <   Ex.Length )
                    {
                        string csvR = string.Format("{0},{1},{2},{3},{4},{5}", Strategy[i], Date[En[i]], Close[En[i]], Date[Ex[i]], Close[Ex[i]], Type[i] * (Close[Ex[i]] - Close[En[i]]));
                        stream.WriteLine(csvR);
                        Console.WriteLine( Strategy[i]+"......."+ Date[En[i]] + "......." + Close[En[i]] + "......." + Date[Ex[i]] + "......." + Close[Ex[i]] + "......." + Type[i] * (Close[Ex[i]] - Close[En[i]]));
                    }
                    else
                    {
                        string csvR = string.Format("{0},{1},{2},{3},{4},{5}","Open"+ Strategy[i], Date[En[i]], Close[En[i]], Date[Date.Length-1], Close[Date.Length - 1], Type[i] * (Close[Date.Length - 1] - Close[En[i]]));
                        stream.WriteLine(csvR);
                        Console.WriteLine("Open" + Strategy[i] + "..." + Date[En[i]] + "......." + Close[En[i]] + "......." + Date[Date.Length - 1] + "......." + Close[Date.Length - 1] + "......." + Type[i] * (Close[Date.Length - 1] - Close[En[i]]));
                    }               
                }
            }
        }
    }
}
