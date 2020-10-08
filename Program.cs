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
    class Program
    {
        static void Main(string[] args) 
        {
            var reader = new StreamReader(File.OpenRead(@"C:\Users\ragin.pm\Documents\Gold.csv"));
            List<string> listA = new List<string>();
            List<string> listB = new List<string>();
            List<string> listC = new List<string>();
            float[] asd = new float[10];
            while (!reader.EndOfStream)
            {

                var line = reader.ReadLine();
                var values = line.Split(',');

                listA.Add(values[0]);
                listB.Add(values[1]);
                listC.Add(values[2]);
            }
            string[] str = listC.ToArray();
            string[] DT = listB.ToArray();
            string[] TR = listA.ToArray();

            float[] Close = Array.ConvertAll(str, float.Parse);
            var movingAverage1 = new MovingAverage();             
            float[] ma1 = movingAverage1.Calc(Close,5);
            float[] ma2 = movingAverage1.Calc(Close,8);      
            float[] ma3 = movingAverage1.Calc(Close,15);
            float[] ma4 = movingAverage1.Calc(Close,30);

            var cr = new CrossOver();
            int[] buy1 = cr.Cross(ma1, ma2);
            int[] sell1 = cr.Cross(ma2, ma1);
            int[] short1 = cr.Cross(ma2, ma1);
            int[] cover1 = cr.Cross(ma1, ma2);

            int[] buy2 = cr.Cross(ma3, ma4);    
               
            int[] sell2 = cr.Cross(ma4, ma3);
            int[] short2 = cr.Cross(ma4, ma3);
            int[] cover2 = cr.Cross(ma3, ma4);

            int[] position = new int[Close.Length];
            Array.Clear(position, 0, position.Length);

            int[] position1 = new int[Close.Length];
            Array.Clear(position1, 0, position.Length);

            //////////////////////////////////////////////
            var Buy1 = new Dictionary<String, Buy>();
            var Short1 = new Dictionary<String, Short>();
            var Sell1 = new Dictionary<String, Sell>();
            var Cover1 = new Dictionary<String, Cover>();

            var Buy2 = new Dictionary<String, Buy>();
            var Short2 = new Dictionary<String, Short>();
            var Sell2 = new Dictionary<String, Sell>();
            var Cover2 = new Dictionary<String, Cover>();
            DateTime Dtp= Convert.ToDateTime(DT[0]); ;
            DateTime Dt;
            int i = 1;
            Entry.Count = 0;
            while(i<Close.Length)
            {
                Dtp = Dtp.AddDays(1);
                Dt = Convert.ToDateTime(DT[i]);
                //Console.WriteLine(Dtp + "-----" + Dt);
                if (Dtp == Dt)
                {
                    //////////////////////////STRATEGY 1///////////////////////////////
                    ///////////////////////////////////////////////////////////////////  
                    if (position[i - 1] == 0 && buy1[i] == 1) { Buy1.Add("Buy1" + i, new Buy(1) { Type = "Buy1", Candle = i, TradeNum = Entry.Count - 1 }); position[i] = 1; }
                    else if (position[i - 1] == 0 && short1[i] == 1) { Short1.Add("Short1" + i, new Short(-1) { Type = "Short1", Candle = i, TradeNum = Entry.Count - 1 }); position[i] = -1; }
                    else if (position[i - 1] == 1 && sell1[i] == 1)
                    {
                        position[i] = -1;
                        Sell1.Add("Sell1" + i, new Sell(Buy1.Values.Last().TradeNum) { Type = "Sell1", Candle = i, TradeNum = Buy1.Values.Last().TradeNum });
                        Short1.Add("Short1" + i, new Short(-1) { Type = "Short1", Candle = i, TradeNum = Entry.Count - 1 });
                    }
                    else if (position[i - 1] == -1 && cover1[i] == 1)
                    {
                        position[i] = 1;
                        Cover1.Add("Cover1" + i, new Cover(Short1.Values.Last().TradeNum) { Type = "Cover1", Candle = i, TradeNum = Short1.Values.Last().TradeNum });
                        Buy1.Add("Buy1" + i, new Buy(1) { Type = "Buy1", Candle = i, TradeNum = Entry.Count - 1 });
                    }
                    else { position[i] = position[i - 1]; }

                    /////////////////////////STRATEGY 2//////////////////////////////
                    /////////////////////////////////////////////////////////////////
                    if (position1[i - 1] == 0 && buy2[i] == 1) { Buy2.Add("Buy1" + i, new Buy(1) { Type = "Buy2", Candle = i, TradeNum = Entry.Count - 1 }); position1[i] = 1; }
                    else if (position1[i - 1] == 0 && short2[i] == 1) { Short2.Add("Short1" + i, new Short(1) { Type = "Short2", Candle = i, TradeNum = Entry.Count - 1 }); position1[i] = -1; }
                    else if (position1[i - 1] == 1 && sell2[i] == 1)
                    {
                        position1[i] = -1;
                        Sell2.Add("Sell1" + i, new Sell(Buy2.Values.Last().TradeNum) { Type = "Sell2", Candle = i, TradeNum = Buy2.Values.Last().TradeNum });
                        Short2.Add("Short1" + i, new Short(-1) { Type = "Short2", Candle = i, TradeNum = Entry.Count - 1 });
                    }
                    else if (position1[i - 1] == -1 && cover2[i] == 1)
                    {
                        position1[i] = 1;
                        Cover2.Add("Cover1" + i, new Cover(Short2.Values.Last().TradeNum) { Type = "Cover2", Candle = i, TradeNum = Short2.Values.Last().TradeNum });
                        Buy2.Add("Buy1" + i, new Buy(1) { Type = "Buy2", Candle = i, TradeNum = Entry.Count - 1 });
                    }
                    else { position1[i] = position1[i - 1]; }
                    ////////
                    ///
                    i = i + 1;
                }
            }
            var Rep = new CalcTrade();
            Rep.Report(Close, DT);                      
        }
    }
}
