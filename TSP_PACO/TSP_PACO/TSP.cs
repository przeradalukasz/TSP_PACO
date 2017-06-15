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
        public static List<Town> Towns { get; set; }

        private static readonly Random rnd = new Random(DateTime.Now.Millisecond);

        public static readonly double Alpha = 1; 
        public static readonly double Beta = 3; 
        public static readonly double Q0 = 0.9; 
        public static readonly double Rho = 0.1; 
        public static readonly double Ksi = 0.1; 
        public static readonly int NumberOfAnts = 300;
        public static readonly int Iterations = 100;

        static void Main()
        {
            Towns = Utils.LoadTownsDataFromJson(@"C:\MagisterkaDane\Dane\Sizes\100\average.json");

            AdjacencyMatrix = Utils.LoadCrispDistanceDataFromJson(@"C:\MagisterkaDane\Dane\crispAdjacencyMatrixNewHampshire.json");
            AdjacencyMatrix = Utils.CutUnnecesssaryElements(AdjacencyMatrix, Towns);
            Towns = Utils.OrderTowns(Towns);

            PACO paco = new PACO(AdjacencyMatrix,Towns.ToArray(),Alpha,Beta,Q0,Rho,Ksi,rnd,NumberOfAnts, Iterations);
            paco.Start();
            Ant ant = paco.GlobalBestDistance;
            Ant ant2 = paco.GlobalBestUnbalancing;

            List<Ant> allTimeBest = paco.CollectionOfFirstFronts;
            List<Ant> allTimeBestFirstFront = paco.GetNondominatedIndividuals(paco.CollectionOfFirstFronts.ToArray());


            Utils.SaveResultsToCsv2(@"C:\MagisterkaDane\Wyniki\105allTimeBestAnt300.csv", allTimeBest.ToArray());
            Utils.SaveResultsToCsv2(@"C:\MagisterkaDane\Wyniki\105allTimeBestFrontAnt300.csv", allTimeBestFirstFront.ToArray());

            //while (population.Generations < 50)
            //{
            //    population.NaturalSelection();
            //    population.Generate();
            //}
            //population.NaturalSelection();

            //Utils.SaveResultsToCsv(@"C:\Users\ronal_000\Desktop\plox.csv", population);


            //Path best = population.Fronts[0].OrderBy(x => x.Distance).First();


            //Utils.DrawResultPath(@"C:\Users\ronal_000\Desktop\dj2.bmp", best.Towns);
            Console.WriteLine("Press any key...");
            Console.ReadLine();

        }
    }
}
