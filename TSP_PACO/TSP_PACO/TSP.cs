using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TSP_PACO
{
    public class TSP
    {
        public static double[,] AdjacencyMatrix { get; set; }
        public static Town[] Towns { get; set; }

        private static readonly Random rnd = new Random(DateTime.Now.Millisecond);

        public static readonly double Alpha = 1; 
        public static readonly double Beta = 3; 
        public static readonly double Q0 = 0.9; 
        public static readonly double Rho = 0.1; 
        public static readonly double Ksi = 0.1; 
        public static readonly int NumberOfAnts = 10;

        static void Main()
        {
            Towns = Utils.LoadTownsData(@"C:\dj.csv");

            AdjacencyMatrix = Utils.LoadDistanceData(Towns);

            PACO paco = new PACO(AdjacencyMatrix,Towns,Alpha,Beta,Q0,Rho,Ksi,rnd,NumberOfAnts,10);
            paco.Start();
            List<Town> best = paco.GlobalBestPath;
            double lol = paco.NearestNeighbour();

            //while (population.Generations < 50)
            //{
            //    population.NaturalSelection();
            //    population.Generate();
            //}
            //population.NaturalSelection();

            //Utils.SaveResultsToCsv(@"C:\Users\ronal_000\Desktop\plox.csv", population);


            //Path best = population.Fronts[0].OrderBy(x => x.Distance).First();


            //Utils.DrawResultPath(@"C:\Users\ronal_000\Desktop\dj2.bmp", best.Towns);

            Console.ReadLine();

        }
    }
}
