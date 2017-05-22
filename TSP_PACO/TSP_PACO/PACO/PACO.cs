using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TSP_PACO
{
    class PACO
    {
        private readonly Random _rnd;
        private readonly double _alpha; //heuristic parameter
        private readonly double _beta; //heuristic parameter
        private readonly double _q0; //control parameter for random proportional
        private readonly double _rho; //evaporation coefficient
        private readonly double _ksi; //local_pheromone
        private readonly int _numberOfAnts;
        private readonly int _iterations;
        private readonly int Q = 10; //constant representing the amount of pheromone an ant put on the path
        public Town[] Towns;
        private readonly int numberOfCities;
        public double[,] AdjacencyMatrix;
        public double[,] PheromoneMatrix;
        public List<Town> GlobalBestPath;
        public double GlobalBestPathLength;


        public Ant[] Ants;

        public PACO(double[,] adjacencyMatrix, Town[] towns, double alpha, double beta, double q0, double rho, double ksi, Random rnd, int numberOfAnts, int iterations)
        {
            _alpha = alpha;
            _beta = beta;
            _q0 = q0;
            _rho = rho;
            _ksi = ksi;
            _rnd = rnd;
            _numberOfAnts = numberOfAnts;
            _iterations = iterations;
            Towns = towns;
            AdjacencyMatrix = adjacencyMatrix;
            numberOfCities = towns.Length;
            PheromoneMatrix = new double[numberOfCities + 1, numberOfCities + 1];
            Ants = new Ant[numberOfAnts];
            GlobalBestPath = new List<Town>();
            GlobalBestPathLength = 0;
        }

        public double CalcDistance(List<Town> tour)
        {
            double result = 0.0;
            for (int i = 0; i < tour.Count - 1; i++)
            {
                result += AdjacencyMatrix[tour[i].Id, tour[i + 1].Id];
            }
            return result;
        }

        //nearest neighbor for initialization pheromone
        public double NearestNeighbour()
        {
            List<Town> path = new List<Town>();
            List<Town> towns = Towns.ToList();
            int remove = 0;
            double tourLength = 0;

            Town nextTown = null;
            Town startingCity = towns.Last();
            path.Add(startingCity);
            towns.Remove(startingCity);

            while (towns.Count > 0)
            {
                double minDistance = AdjacencyMatrix[startingCity.Id, towns[0].Id];
                remove = 0;
                for (int i = 0; i < towns.Count; i++)
                {
                    double distance = AdjacencyMatrix[startingCity.Id, towns[i].Id];
                    if (Math.Abs(distance) > 0.00000001 && distance < minDistance)
                    {
                        minDistance = distance;
                        nextTown = towns[i];
                        remove = i;
                    }
                }
                startingCity = nextTown;
                towns.RemoveAt(remove);
                path.Add(nextTown);
                tourLength = tourLength + minDistance;
            }

            path.Add(path[0]);
            tourLength = tourLength + AdjacencyMatrix[nextTown.Id, path[0].Id];

            return tourLength;
        }

        public void InitializePheromone(double tau0)
        {
            for (int i = 0; i < numberOfCities; i++)
            {
                for (int j = 0; j < numberOfCities; j++)
                {
                    if (i == j)
                    {
                        PheromoneMatrix[i, j] = 0;
                    }
                    else
                    {
                        PheromoneMatrix[i, j] = tau0;
                    }

                }

            }
        }

        public double CalculateTauEtha(Town currentTown, out List<double> tauEtha, List<Town> citiesLeft)
        {
            tauEtha = new List<double>();
            double total = 0;
            for (int i = 0; i < citiesLeft.Count; i++)
            {
                double tauEthaVal;
                try
                {
                    tauEthaVal = Math.Pow(PheromoneMatrix[currentTown.Id, i], _alpha) *
                                 Math.Pow(1.0 / AdjacencyMatrix[currentTown.Id, citiesLeft[i].Id], _beta);
                    if (double.IsPositiveInfinity(tauEthaVal))
                    {
                        tauEthaVal = 0;
                    }
                }
                catch (DivideByZeroException e)
                {
                    tauEthaVal = 0;
                    throw;
                }
                tauEtha.Add(tauEthaVal);
                total = total + tauEthaVal;
            }
            return total;
        }

        public Town FindNextCity(Town currentTown, List<Town> citiesLeft)
        {
            List<double> tauEtha;

            double totalTauEtha = CalculateTauEtha(currentTown, out tauEtha, citiesLeft);

            if (_rnd.NextDouble() < _q0)
            {
                double argmax = tauEtha.Max();
                return citiesLeft[tauEtha.IndexOf(argmax)];
            }
            else
            {
                double roulette = 0;
                for (int i = 0; i < citiesLeft.Count; i++)
                {
                    roulette = roulette + tauEtha[i];
                    if (_rnd.Next((int)totalTauEtha) < roulette)
                        return citiesLeft[i];
                }
            }
            return null;
        }

        public void TourConstruction(Ant ant)
        {
            List<Town> towns = Towns.ToList();
            Town curretnTown = ant.Tour[0];
            towns.Remove(curretnTown);

            for (int i = 0; i < numberOfCities - 2; i++)
            {
                Town nextTown = FindNextCity(curretnTown, towns);
                ant.Tour.Add(nextTown);
                curretnTown = nextTown;
                towns.Remove(nextTown);
            }
            ant.Tour.Add(towns.Last());
            ant.Tour.Add(ant.Tour[0]);
            ant.TourLength = CalcDistance(ant.Tour);
        }

        public void LocalPheromoneUpdate(Ant ant, double tau0)
        {
            for (int i = 0; i < ant.Tour.Count - 1; i++)
            {
                int current = ant.Tour[i].Id;
                int next = ant.Tour[i + 1].Id;

                PheromoneMatrix[current, next] = ((1 - _ksi) * PheromoneMatrix[current, next]) + (_ksi * tau0);
                PheromoneMatrix[next, current] = PheromoneMatrix[current, next];
            }

        }

        public void GlobalPheromoneUpdate(List<Town> globalBestTour, double globalBestTourLenght)
        {
            for (int i = 0; i < globalBestTour.Count - 1; i++)
            {
                int current = globalBestTour[i].Id;
                int next = globalBestTour[i + 1].Id;

                PheromoneMatrix[current, next] = ((1 - _rho) * PheromoneMatrix[current, next]) + (Q * (1 / globalBestTourLenght));
                PheromoneMatrix[next, current] = PheromoneMatrix[current, next];
            }
        }

        public Tuple<List<Town>, double> LocalSearch(int andId, List<Town> antTour, double antTourLength)
        {
            while (true)
            {
                double best = antTourLength;

                for (int i = 0; i < antTour.Count - 1; i++)
                {
                    for (int j = i + 1; j < antTour.Count; j++)
                    {
                        List<Town> newAntTour = new List<Town>(antTour);
                        int k = i;
                        int l = j;

                        while (k < l)
                        {
                            Town swap = newAntTour[k];
                            newAntTour[k] = newAntTour[l];
                            newAntTour[l] = swap;

                            if (k == 0)
                            {
                                newAntTour[antTour.Count - 1] = newAntTour[k];
                            }
                            if (l == antTour.Count-1)
                            {
                                newAntTour[0] = newAntTour[l];
                            }

                            k = k + 1;
                            l = l - 1;
                        }

                        double newAntTourLength = CalcDistance(newAntTour);

                        if (newAntTourLength < antTourLength)
                        {
                            antTourLength = newAntTourLength;
                            antTour = newAntTour;
                        }
                    }
                }
                if (Math.Abs(best - antTourLength) < 0.000000001)
                {
                    return new Tuple<List<Town>, double>(antTour, antTourLength);
                }
            }
        }

        public void InitializeTours(List<Town> bestTour)
        {
            bestTour = new List<Town>();
            for (int i = 0; i < _numberOfAnts; i++)
            {
                Ants[i].Tour = new List<Town>();
                Ants[i].TourLength = 0;
                Ants[i].Tour.Add(Towns[_rnd.Next(numberOfCities)]);
            }
        }

        public void SystemStart(double tau0)
        {
            InitializePheromone(tau0);
            List<Town> bestTour = new List<Town>();

            for (int i = 0; i < _iterations; i++)
            {
                double bestTourLength = 0;
                InitializeTours(bestTour);
                for (int j = 0; j < _numberOfAnts; j++)
                {
                    TourConstruction(Ants[j]);
                    LocalPheromoneUpdate(Ants[j], tau0);

                    var localSearchResult = LocalSearch(Ants[j].Id, Ants[j].Tour, Ants[j].TourLength);
                    Ants[j].Tour = localSearchResult.Item1;
                    Ants[j].TourLength = localSearchResult.Item2;

                    if (Math.Abs(bestTourLength) < 0.0000000001 || bestTourLength > Ants[j].TourLength)
                    {
                        bestTourLength = Ants[j].TourLength;
                        bestTour = Ants[j].Tour;
                    }
                }
                if (Math.Abs(GlobalBestPathLength) < 0.0000000001 || GlobalBestPathLength > bestTourLength)
                {
                    GlobalBestPathLength = bestTourLength;
                    GlobalBestPath = bestTour;
                }
                GlobalPheromoneUpdate(GlobalBestPath, GlobalBestPathLength);
            }
        }

        public void Start()
        {
            double tau0 = 1.0 / numberOfCities * NearestNeighbour();

            for (int i = 0; i < _numberOfAnts; i++)
            {
                Ants[i] = new Ant(i);
            }
            SystemStart(tau0);
        }


    }
}
