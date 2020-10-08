using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ConsoleApp42
{
    public abstract class Trade
    {
        public string Type { get; set; }
        public int Candle { get; set; }
        public int TradeNum { get; set; }
        public int globalId { get; set; }
        public int SomeValue = 0;
        public static List<Trade> EntryList = new List<Trade>();
        public static List<Trade> ExitList = new List<Trade>();

    }
    public class Entry : Trade
    {
        public static int Count { get; set; }
        public Entry()
        {
            globalId = Count; Count++;
            EntryList.Add(this);
        }
    }
    public class Exit : Trade
    {
        public Exit()
        {

        }
    }

    public class Buy : Entry
    {
        public Buy(int i) : base()
        {
            this.SomeValue = i;
        }
    }
    public class Short : Entry
    {
        public Short(int i) : base()
        {
            this.SomeValue = i;
        }
    }
    public class Sell : Exit
    {
        public Sell(int i) : base()
        {
            this.SomeValue = i;
            globalId = i;
            if (i - (ExitList.Count() - 1) > 1 || i - ExitList.Count() == 0)
            {
                ExitList.Add(this);
            }
            else
            {
                ExitList.Insert(i, this);
            }
        }
    }

    public class Cover : Exit
    {
        public Cover(int i) : base()
        {
            this.SomeValue = i;
            globalId = i;
            if (i - (ExitList.Count() - 1) > 1 || i - ExitList.Count() == 0)
            {
                ExitList.Add(this);
            }
            else
            {
                ExitList.Insert(i, this);
            }
        }
    }
}
