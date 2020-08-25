using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using MathNet;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra;

namespace GARCH
{
    public class funct
    {
        public int Roll;
        public double[] MA(double[] Close) //Function For Moving Average
        {
            double[] Ma = new double[Close.Length];
            for (int i = Roll - 1; i < Close.Length; i++)
            {
                double sum = 0;

                for (int j = 0; j < Roll; j++)
                {
                    sum = sum + Close[i - j];
                }
                Ma[i] = sum / Roll;
            }
            return Ma;
        }
        public double[] Var(double[] eps)//Function For Varience
        {
            double[] MA1 = MA(eps);
            double[] var = new double[eps.Length];
            for (int i = Roll - 1; i < eps.Length; i++)
            {
                double sum = 0;
                for (int j = 0; j < Roll; j++)
                {
                    sum = sum + (MA1[i] - eps[i - j]) * (MA1[i] - eps[i - j]);
                }
                var[i] = sum / Roll;
            }
            return var;
        }
        public double[] autocorr(double[] ret)//Function For Aurocorrelation
        {
            double[] corr = new double[ret.Length];
            double[] MA1 = MA(ret);
            for (int m = Roll + 1; m < ret.Length; m++)
            {
                double multi = 0;
                for (int j = m - Roll + 1 + 1; j < m + 1; j++)
                {
                    multi = (ret[j] - MA1[m]) * (ret[j - 1] - MA1[m]) + multi;
                }
                double multi2 = 0;
                for (int j = m - Roll + 1; j < m + 1; j++)
                {
                    multi2 = (ret[j] - MA1[m]) * (ret[j] - MA1[m]) + multi2;
                }
                corr[m] = multi / multi2;
            }
            return corr;
        }
        public double Fn(int a)
        {
            int R = 0;
            if (a <= 0) { R = 0; }
            else { R = 1; }
            return R;
        }
    }
    //////////////////////////////////////////////////////////
    ///////////////////////GARCH Calculation//////////////////
    //////////////////////////////////////////////////////////
    public class GARCH
    {
        public int Roll;
        public double alpha;
        public int r;
        public int q;
        public int p;
        public double[,] Calc(double[] Close)
        {

            double[] ret = new double[Close.Length];
            for (int i = 1; i < Close.Length; i++)
            {
                ret[i] = Close[i] - Close[i - 1];
            }
            double[,] b = new double[Close.Length, r + 1];
            double[,] bi = new double[Close.Length, r + 1];

            double[,] x = new double[Close.Length, r + 1];

            double[,] z = new double[Close.Length, q + p];
            double[,] zi = new double[Close.Length, q + p];

            double[,] w = new double[Close.Length, q + p];
            double[,] wi = new double[Close.Length, q + p];

            double[,] F_theta = new double[Close.Length, 1];
            double[] theta = new double[r + 1 + q + p];
            double[] s = new double[r + 1 + q + p];
            double[] thetai = new double[r + 1 + q + p];
            double[] L_theta = new double[r + 1 + q + p];
            double[] L_theta1 = new double[r + 1 + q + p];
            //double[] L_theta2 = new double[r + 1 + q + 1 + p];
            //double[,] H       = new double[r + 1 + q + 1 + p, r + 1 + q + 1 + p];
            Vector<double> L_theta2 = Vector<double>.Build.Dense(r + 1 + q + p);
            Matrix<double> H = Matrix<double>.Build.Dense(r + 1 + q + p, r + 1 + q + p);
            double[] eps = new double[Close.Length];
            double[] epsi = new double[Close.Length];

            double[] h = new double[Close.Length];
            double[] hi = new double[Close.Length];

            double[] ht = new double[Close.Length];
            double[] hti = new double[Close.Length];

            double[] b0 = new double[Close.Length];
            double[] b0i = new double[Close.Length];

            var Fn = new funct();
            Fn.Roll = Roll;
            double[] b1 = Fn.autocorr(ret); //AR(1) Coefficient
            double[] b1i = Fn.autocorr(ret);
            double[] mu = Fn.MA(ret);

            double[,] D_h_w = new double[Close.Length, q + p];
            double[,] D_h_b = new double[Close.Length, r + 1];
            double[,] D_l_w = new double[Close.Length, q + p];
            double[,] D_l_b = new double[Close.Length, r + 1];

            double[,] outp = new double[Close.Length, 7];
            double[] lk = new double[2]; //Likelihoosd Function For Convergence

            double[] sum1 = new double[r + 1];
            double sum = 0;
            double sum2 = 0;
            int val = 0;

            for (int i = 1; i < Close.Length; i++)
            {
                b0[i] = mu[i] * (1 - b1[i]);
                eps[i] = ret[i] - (b0[i] + b1[i] * ret[i - 1]);
                epsi[i] = ret[i] - (b0[i] + b1[i] * ret[i - 1]);

                b[i, 0] = b0[i];
                b[i, 1] = b1[i];

                bi[i, 0] = b0[i];
                bi[i, 1] = b1[i];

                ht[i] = eps[i] * eps[i];
                hti[i] = eps[i] * eps[i];

                outp[i, 0] = ht[i];

                x[i, 0] = 1;
                x[i, 1] = ret[i - 1];
            }
            h = Fn.Var(eps);
            hi = Fn.Var(epsi);
            for (int i = 1; i < Close.Length; i++)
            {
                //z[i, 0] = 1;
                z[i, 0] = ht[i - 1];
                z[i, 1] = h[i - 1];

                //zi[i, 0] = 1;
                zi[i, 0] = ht[i - 1];
                zi[i, 1] = h[i - 1];

                outp[i, 1] = h[i];
            }
            ///////////////////////START/////////////////////
            /////////////////////////////////////////////////
            for (int t = 3 * Roll - 2; t < Close.Length; t++)
            {
                /////////////////Initialising/////////////////////
                for (int i = t - Roll - p - q + 1; i < t + 1; i++)
                {
                    z[i, 0] = zi[i, 0];
                    z[i, 1] = zi[i, 1];


                    eps[i] = epsi[i];
                    ht[i] = hti[i];
                    h[i] = hi[i];

                    b[i, 0] = bi[i, 0];
                    b[i, 1] = bi[i, 1];
                }
                lk[0] = 0;
                lk[1] = 0;
                sum2 = 0;
                for (int j = 1; j < r + 1; j++)
                {
                    sum1[j] = 0;
                }
                for (int j = 0; j < r + 1 + q + p; j++)
                {
                    if (j < r + 1)
                    {
                        theta[j] = bi[t, j];
                    }
                    else
                    {
                        theta[j] = 0; ;
                        w[t, j - (r + 1)] = 0; ;
                    }

                }

                for (int j = 0; j < r + 1 + q + p; j++)
                {
                    L_theta1[j] = 0;
                    L_theta[j] = 0;
                    L_theta2[j] = 0;
                }
                for (int j = 0; j < r + 1 + q + p; j++)
                {
                    for (int k = 0; k < r + 1 + q + p; k++)
                    {
                        H[j, k] = 0;
                    }
                }

                for (int i = t - Roll - p - q + 1; i < t + 1; i++)
                {
                    for (int j = 0; j < q + p; j++)
                    {
                        if (j < r + 1)
                        {
                            D_h_b[i, j] = 0;
                            D_l_b[i, j] = 0;
                        }
                        else
                        {
                            D_h_w[i, j - r - 1] = 0;
                            D_l_w[i, j - r - 1] = 0;
                        }
                    }
                }

                for (int j = 0; j < r + 1; j++)
                {
                    sum1[j] = 0;
                }
                for (int i = t - Roll + 1; i < t + 1; i++)
                {
                    for (int j = 0; j < r + 1; j++)
                    {
                        if (i < t)
                        {
                            sum1[j] = sum1[j] + eps[i] * x[i, j];
                        }
                        else
                        {
                            sum1[j] = (1 / Roll) * sum1[j] + eps[i] * x[i, j];
                        }
                    }
                }
                //////////Start Roll Initialization///////
                //////////////////////////////////////////
                for (int i = t - Roll + 1; i < t + 1; i++)
                {
                    for (int j = 0; j < q + p; j++)
                    {
                        if (i == t - Roll + 1)
                        {
                            D_h_w[i, j] = 0;
                        }
                        else
                        {
                            sum = 0;
                            for (int k = 1; k < p + 1; k++)
                            {
                                sum = sum + w[t, k] * D_h_w[i - k, j];
                            }
                            D_h_w[i, j] = z[i, j] + sum;
                        }
                        D_l_w[i, j] = 0.5d * (1 / h[i]) * D_h_w[i, j] * ((ht[i] / h[i]) - 1);
                        //Console.WriteLine(z[i,j]);
                    }

                    for (int j = 0; j < r + 1; j++)
                    {
                        sum = 0;
                        for (int k = 1; k < q + 1; k++)
                        {
                            sum = sum + w[t, k] * Math.Pow(sum1[j], 1 - (Fn.Fn(val - k))) * Math.Pow((x[i - k, j] * eps[i - k]), (Fn.Fn(val - k)));
                        }
                        if (i == t - Roll + 1)
                        {
                            D_h_b[i, j] = -2 * sum1[j];
                        }
                        else
                        {
                            sum2 = 0;
                            for (int k = 1; k < p + 1; k++)
                            {
                                sum2 = sum2 + w[t, k] * D_h_b[i - k, j];
                            }
                            D_h_b[i, j] = sum + sum2;
                        }
                        D_l_b[i, j] = ((eps[i] * x[i, j]) / h[i]) + (0.5d) * (1 / h[i]) * D_h_b[i, j] * ((ht[i] / h[i]) - 1);
                    }

                    for (int j = 0; j < r + 1 + q + p; j++)
                    {
                        if (j < r + 1)
                        {
                            L_theta[j] = D_l_b[i, j] + L_theta[j];

                            L_theta1[j] = D_l_b[i, j];
                            L_theta2[j] = L_theta[j];
                        }
                        else
                        {
                            L_theta[j] = D_l_w[i, j - (r + 1)] + L_theta[j];

                            L_theta1[j] = D_l_w[i, j - (r + 1)];
                            L_theta2[j] = L_theta[j];
                        }
                    }
                    for (int j = 0; j < r + 1 + q + p; j++)
                    {
                        for (int k = 0; k < r + 1 + q + p; k++)
                        {
                            H[j, k] = L_theta1[j] * L_theta1[k] + H[j, k];
                        }
                    }
                    val = val + 1;
                }
                H = H.Inverse();
                //Console.Write(L_theta2);
                double er = 10;
                ///////Itration///////
                //////////////////////
                while (er > 0.001)//for(int itr = 0; itr<10; itr++)
                {
                    ///////////////////////////////////////////////////////
                    H.Multiply(L_theta2, L_theta2);
                    //Console.Write(L_theta2);
                    for (int j = 0; j < r + 1 + q + p; j++)
                    {
                        thetai[j] = theta[j];
                        theta[j] = theta[j] + alpha * (L_theta2[j]);//Itration algorith//
                    }

                    for (int j = 0; j < r + 1 + q + p; j++)
                    {
                        if (j < r + 1)
                        {
                            b[t, j] = theta[j];
                        }
                        else
                        {
                            w[t, j - (r + 1)] = theta[j];
                        }
                    }

                    for (int i = t - Roll - p + 1; i < t + 1; i++)
                    {
                        eps[i] = ret[i] - (b[t, 0] + b[t, 1] * ret[i - 1]);
                        ht[i] = eps[i] * eps[i];
                    }

                    for (int m = 0; m < p + 1; m++)
                    {
                        sum = 0;
                        for (int i = t - Roll - m + 1; i < t + 1 - m; i++)
                        {
                            sum = sum + eps[i];
                        }
                        sum2 = 0;
                        sum2 = sum / Roll;
                        sum = 0;
                        for (int i = t - Roll - m + 1; i < t + 1 - m; i++)
                        {
                            sum = sum + (eps[i] - sum2) * (eps[i] - sum2);
                        }
                        h[t - m] = sum / Roll;
                    }
                    outp[t, 4] = h[t] * (1 - (w[t, 0] + w[t, 1]));
                    z[t, 0] = ht[t - 1];
                    z[t, 1] = h[t - 1];

                    sum = 0;
                    for (int i = t - Roll + 1; i < t + 1; i++)
                    {
                        sum = sum - 1 * ((float)Math.Log(h[t]) + (ht[t]) / h[t]);
                    }
                    lk[1] = 0;
                    lk[1] = sum;
                    er = Math.Abs(lk[1] - lk[0]);
                    lk[0] = lk[1];
                    sum = 0;
                    ///////////////////////////////////////////////////////
                    for (int j = 0; j < r + 1 + q + p; j++)
                    {
                        L_theta1[j] = 0;
                        L_theta[j] = 0;
                        L_theta2[j] = 0;
                    }
                    for (int i = t - Roll - p - q + 1; i < t + 1; i++)
                    {
                        for (int j = 0; j < q + p; j++)
                        {
                            if (j < r + 1)
                            {
                                D_h_b[i, j] = 0;
                                D_l_b[i, j] = 0;
                            }
                            else
                            {
                                D_h_w[i, j - r - 1] = 0;
                                D_l_w[i, j - r - 1] = 0;
                            }
                        }
                    }

                    for (int j = 0; j < r + 1; j++)
                    {
                        sum1[j] = 0;
                    }
                    for (int i = t - Roll + 1; i < t + 1; i++)
                    {
                        for (int j = 0; j < r + 1; j++)
                        {
                            if (i < t)
                            {
                                sum1[j] = sum1[j] + eps[i] * x[i, j];
                            }
                            else
                            {
                                sum1[j] = (1 / Roll) * sum1[j] + eps[i] * x[i, j];
                            }
                        }
                    }
                    ///////////////////Start Roll ////////////
                    //////////////////////////////////////////
                    for (int i = t - Roll + 1; i < t + 1; i++)
                    {
                        for (int j = 0; j < q + p; j++)
                        {
                            if (i == t - Roll + 1)
                            {
                                D_h_w[i, j] = 0;
                            }
                            else
                            {
                                sum = 0;
                                for (int k = 1; k < p + 1; k++)
                                {
                                    sum = sum + w[t, k] * D_h_w[i - k, j];
                                }
                                D_h_w[i, j] = z[i, j] + sum;
                            }
                            D_l_w[i, j] = 0.5d * (1 / h[i]) * D_h_w[i, j] * ((ht[i] / h[i]) - 1);
                        }

                        for (int j = 0; j < r + 1; j++)
                        {
                            sum = 0;
                            for (int k = 1; k < q + 1; k++)
                            {
                                sum = sum + w[t, k] * Math.Pow(sum1[j], 1 - (Fn.Fn(val - k))) * Math.Pow((x[i - k, j] * eps[i - k]), (Fn.Fn(val - k)));
                            }
                            if (i == t - Roll + 1)
                            {
                                D_h_b[i, j] = -2 * sum1[j];
                            }
                            else
                            {
                                sum2 = 0;
                                for (int k = 1; k < p + 1; k++)
                                {
                                    sum2 = sum2 + w[t, k] * D_h_b[i - k, j];
                                }
                                D_h_b[i, j] = sum + sum2;
                            }
                            D_l_b[i, j] = ((eps[i] * x[i, j]) / h[i]) + (0.5d) * (1 / h[i]) * D_h_b[i, j] * ((ht[i] / h[i]) - 1);
                        }

                        for (int j = 0; j < r + 1 + q + p; j++)
                        {
                            if (j < r + 1)
                            {
                                L_theta[j] = D_l_b[i, j] + L_theta[j];

                                L_theta1[j] = D_l_b[i, j];
                                L_theta2[j] = L_theta[j];
                            }
                            else
                            {
                                L_theta[j] = D_l_w[i, j - (r + 1)] + L_theta[j];

                                L_theta1[j] = D_l_w[i, j - (r + 1)];
                                L_theta2[j] = L_theta[j];
                            }
                        }
                    }
                    /////////////////////End Roll///////////////

                    //////////////Updation of Hessian///////////////////
                    sum = 0;
                    for (int j = 0; j < r + 1 + q + p; j++)
                    {
                        s[j] = theta[j] - thetai[j];
                        sum = sum + L_theta[j] * s[j];
                        thetai[j] = theta[j];
                    }
                    double rho = 1 / sum;
                    sum = 0;
                    Matrix<double> A = Matrix<double>.Build.Dense(r + 1 + q + p, r + 1 + q + p);
                    Matrix<double> B = Matrix<double>.Build.Dense(r + 1 + q + p, r + 1 + q + p);
                    Matrix<double> I = Matrix<double>.Build.Dense(r + 1 + q + p, r + 1 + q + p);
                    for (int i = 0; i < r + 1 + q + p; i++)
                    {
                        for (int j = 0; j < r + 1 + q + p; j++)
                        {
                            if (i == j) { I[i, j] = 1; }
                            else { I[i, j] = 0; }
                            A[i, j] = I[i, j] - rho * s[i] * L_theta[j];
                            B[i, j] = I[i, j] - rho * L_theta[i] * s[j];
                        }
                    }
                    A.Multiply(H, H);
                    H.Multiply(B, H);
                    //Console.Write(H);
                }
                ////////////End While////////
                outp[t, 2] = b[t, 0];// ht[t];
                outp[t, 3] = b[t, 1];// h[t];

                outp[t, 5] = w[t, 0];
                outp[t, 6] = w[t, 1];
            }
            return outp;
        }
    }
    class Program
    {
        static void Main(string[] args)
        {
            var reader = new StreamReader(File.OpenRead(@"C:\Users\ragin.pm\Documents\Gold.csv"));
            List<string> listA = new List<string>();
            List<string> listB = new List<string>();
            List<string> listC = new List<string>();
            double[] asd = new double[10];
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

            double[] Close = Array.ConvertAll(str, double.Parse);
            var G = new GARCH();
            G.alpha = 1;// 0.01d;
            G.p = 1;
            G.q = 1;
            G.r = 1;
            G.Roll = 30;
            double[,] h = G.Calc(Close);
            var rowCount = h.GetLength(0);
            var colCount = h.GetLength(1);
            double[] var = new double[Close.Length];
            for (int i = 1; i < Close.Length; i++)
            {
                var[i] = h[i, 4] + h[i - 1, 0] * h[i, 5] + h[i - 1, 1] * h[i, 6];
            }
            for (int row = 0; row < rowCount; row++)
            {
                for (int col = 0; col < colCount; col++)
                    Console.Write(String.Format("{0}\t", h[row, col]));
                Console.WriteLine();
            }
            var file = @"C:\Users\ragin.pm\Documents\Gold_1.csv";
            using (var stream = File.CreateText(file))
            {

                for (int i = 0; i < Close.Length; i++)
                {
                    string csvR = string.Format("{0},{1},{2},{3},{4},{5},{6},{7}", h[i, 0].ToString(), h[i, 1].ToString(), h[i, 2].ToString(), h[i, 3].ToString(), h[i, 4].ToString(), h[i, 5].ToString(), h[i, 6].ToString(), var[i].ToString());
                    stream.WriteLine(csvR);
                }
            }
        }
    }
}

