using System;
using System.Collections.Generic;
using System.IO;
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
        public double[,] DistancePheromoneMatrix;
        public double[,] UnbalancingDegreePheromoneMatrix;
        public Ant GlobalBestDistance;
        public Ant GlobalBestUnbalancing;
        public Ant GlobalSecondBestDistance;
        public Ant GlobalSecondBestUnbalancing;
        public List<List<Ant>> Fronts;
        public List<Ant> CollectionOfFirstFronts;


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
            DistancePheromoneMatrix = new double[numberOfCities + 1, numberOfCities + 1];
            UnbalancingDegreePheromoneMatrix = new double[numberOfCities + 1, numberOfCities + 1];
            Ants = new Ant[numberOfAnts];
            GlobalBestDistance = new Ant();
            GlobalBestUnbalancing = new Ant();
            GlobalSecondBestDistance = new Ant();
            GlobalSecondBestUnbalancing = new Ant();
            Fronts = new List<List<Ant>>();
            CollectionOfFirstFronts = new List<Ant>();
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

        public double CalcUnbalancingDegree(List<Town> tour)
        {
            double minDistance = 100000.0;
            double maxDistance = 0.0;

            for (int i = 0; i < tour.Count - 1; i++)
            {
                double currentDistance = AdjacencyMatrix[tour[i].Id, tour[i + 1].Id];
                if (currentDistance < minDistance)
                    minDistance = currentDistance;
                if (currentDistance > maxDistance)
                    maxDistance = currentDistance;
            }
            return maxDistance - minDistance;
        }

        //nearest neighbor for initialization pheromone
        public Tuple<double,double> NearestNeighbourMin()
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
                nextTown = towns[0];
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
            var tourAmplitude = CalcUnbalancingDegree(path);
            return new Tuple<double, double>(tourLength, tourAmplitude);
        }

        public Tuple<double, double> NearestNeighbourMax()
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
                double maxDistance = AdjacencyMatrix[startingCity.Id, towns[0].Id];
                nextTown = towns[0];
                remove = 0;
                for (int i = 0; i < towns.Count; i++)
                {
                    double distance = AdjacencyMatrix[startingCity.Id, towns[i].Id];
                    if (Math.Abs(distance) > 0.00000001 && distance > maxDistance)
                    {
                        maxDistance = distance;
                        nextTown = towns[i];
                        remove = i;
                    }
                }
                startingCity = nextTown;
                towns.RemoveAt(remove);
                path.Add(nextTown);
                tourLength = tourLength + maxDistance;
            }

            path.Add(path[0]);
            tourLength = tourLength + AdjacencyMatrix[nextTown.Id, path[0].Id];
            var tourAmplitude = CalcUnbalancingDegree(path);
            return new Tuple<double, double>(tourLength, tourAmplitude);
        }

        public void InitializePheromone(double distanceTau0, double unbalancingTau0)
        {
            for (int i = 0; i < numberOfCities; i++)
            {
                for (int j = 0; j < numberOfCities; j++)
                {
                    if (i == j)
                    {
                        DistancePheromoneMatrix[i, j] = 0;
                        UnbalancingDegreePheromoneMatrix[i, j] = 0;
                    }
                    else
                    {
                        DistancePheromoneMatrix[i, j] = distanceTau0;
                        UnbalancingDegreePheromoneMatrix[i, j] = unbalancingTau0;
                    }

                }

            }
        }

        public double CalculateTauEtha(Town currentTown, out List<double> tauEtha, List<Town> citiesLeft, List<Town> currentTour)
        {
            
            tauEtha = new List<double>();
            double total = 0;
            for (int i = 0; i < citiesLeft.Count; i++)
            {
                double tauEthaVal;
                try
                {
                    //tauEthaVal = Math.Pow(DistancePheromoneMatrix[currentTown.Id, i], _alpha) *
                    //             Math.Pow(1.0 / AdjacencyMatrix[currentTown.Id, citiesLeft[i].Id], _beta);
                    List<Town> possibleTour = new List<Town>(currentTour);
                    possibleTour.Add(citiesLeft[i]);
                    double futureUnbalanceDegree = CalcUnbalancingDegree(possibleTour);
                    double pheromoneWeight = _rnd.NextDouble(); 
                    tauEthaVal = Math.Pow(pheromoneWeight*DistancePheromoneMatrix[currentTown.Id, i] + (1-pheromoneWeight)*UnbalancingDegreePheromoneMatrix[currentTown.Id, i], _alpha) *
                                 Math.Pow(2.0 / (AdjacencyMatrix[currentTown.Id, citiesLeft[i].Id] + futureUnbalanceDegree), _beta);
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

        public Town FindNextCity(Town currentTown, List<Town> citiesLeft, List<Town> currentTour)
        {
            List<double> tauEtha;

            double totalTauEtha = CalculateTauEtha(currentTown, out tauEtha, citiesLeft, currentTour);
           
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
                    if (_rnd.Next((int)totalTauEtha) <= roulette)
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
                Town nextTown = FindNextCity(curretnTown, towns, ant.Tour);
                ant.Tour.Add(nextTown);
                curretnTown = nextTown;
                towns.Remove(nextTown);
            }
            ant.Tour.Add(towns.Last());
            ant.Tour.Add(ant.Tour[0]);
            ant.TourLength = CalcDistance(ant.Tour);
            ant.UnbalancingDegree = CalcUnbalancingDegree(ant.Tour);
        }

        public void LocalPheromoneUpdate(Ant ant, double distanceTau0, double unbalancingTau0)
        {
            for (int i = 0; i < ant.Tour.Count - 1; i++)
            {
                int current = ant.Tour[i].Id;
                int next = ant.Tour[i + 1].Id;

                DistancePheromoneMatrix[current, next] = ((1 - _ksi) * DistancePheromoneMatrix[current, next]) + (_ksi * distanceTau0);
                DistancePheromoneMatrix[next, current] = DistancePheromoneMatrix[current, next];

                UnbalancingDegreePheromoneMatrix[current, next] = ((1 - _ksi) * UnbalancingDegreePheromoneMatrix[current, next]) + (_ksi * unbalancingTau0);
                UnbalancingDegreePheromoneMatrix[next, current] = UnbalancingDegreePheromoneMatrix[current, next];
            }

        }

        public void GlobalPheromoneUpdate(Ant globalBestDistance, Ant globalBestUnbalancing, Ant globalSecondBestDistance, Ant globalSecondBestUnbalancing)
        {
            for (int i = 0; i < globalBestDistance.Tour.Count - 1; i++)
            {
                int current = globalBestDistance.Tour[i].Id;
                int next = globalBestDistance.Tour[i + 1].Id;

                DistancePheromoneMatrix[current, next] = ((1 - _rho) * DistancePheromoneMatrix[current, next]) + (_rho * (Q/globalBestDistance.TourLength));
                DistancePheromoneMatrix[next, current] = DistancePheromoneMatrix[current, next];
            }

            for (int i = 0; i < globalBestUnbalancing.Tour.Count - 1; i++)
            {
                int current = globalBestUnbalancing.Tour[i].Id;
                int next = globalBestUnbalancing.Tour[i + 1].Id;

                UnbalancingDegreePheromoneMatrix[current, next] = ((1 - _rho) * UnbalancingDegreePheromoneMatrix[current, next]) + (_rho * (Q / globalBestDistance.UnbalancingDegree));
                UnbalancingDegreePheromoneMatrix[next, current] = UnbalancingDegreePheromoneMatrix[current, next];
            }

            for (int i = 0; i < globalSecondBestDistance.Tour.Count - 1; i++)
            {
                int current = globalSecondBestDistance.Tour[i].Id;
                int next = globalSecondBestDistance.Tour[i + 1].Id;

                DistancePheromoneMatrix[current, next] = ((1 - _rho) * DistancePheromoneMatrix[current, next]) + (_rho * 0.5 * (Q / globalBestDistance.TourLength));
                DistancePheromoneMatrix[next, current] = DistancePheromoneMatrix[current, next];
            }

            for (int i = 0; i < globalSecondBestUnbalancing.Tour.Count - 1; i++)
            {
                int current = globalSecondBestUnbalancing.Tour[i].Id;
                int next = globalSecondBestUnbalancing.Tour[i + 1].Id;

                UnbalancingDegreePheromoneMatrix[current, next] = ((1 - _rho) * UnbalancingDegreePheromoneMatrix[current, next]) + (_rho * 0.5 *(Q / globalBestDistance.UnbalancingDegree));
                UnbalancingDegreePheromoneMatrix[next, current] = UnbalancingDegreePheromoneMatrix[current, next];
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

        public void InitializeTours()
        {
            for (int i = 0; i < _numberOfAnts; i++)
            {
                Ants[i].Tour = new List<Town>();
                Ants[i].TourLength = 0;
                Ants[i].UnbalancingDegree = 0;
                Ants[i].Tour.Add(Towns[_rnd.Next(numberOfCities)]);
            }
        }

        public Tuple<double, double> NormalizeValues()
        {
            var nearestMax = NearestNeighbourMax();
            var nearestMin = NearestNeighbourMin();
            var distanceNormalizedValue = (nearestMin.Item1 - 0) / (nearestMax.Item1 - 0);
            var unbalancingNormalizedValue = (nearestMin.Item2 - 0) / (nearestMax.Item2 - 0);
            return new Tuple<double, double>(distanceNormalizedValue,unbalancingNormalizedValue);
        }

        public List<Ant> GetNondominatedIndividuals(Ant[] population)
        {
            var pNondominated = new List<Ant>();
            foreach (var individualToCheck in population)
            {
                bool isDominated = population.Where(individual => !individualToCheck.Equals(individual)).Any(individual => Dominates(individualToCheck, individual));

                if (!isDominated)
                {
                    pNondominated.Add(individualToCheck);
                }
            }
            return pNondominated;
        }

        public bool Dominates(Ant individualToCheck, Ant individual)
        {
            bool betterForAllCriteriums = true;

            if (individualToCheck.TourLength < individual.TourLength || individualToCheck.UnbalancingDegree < individual.UnbalancingDegree)
            {
                betterForAllCriteriums = false;
            }

            return betterForAllCriteriums;
        }

        public void GetFronts(Ant[] ants)
        {
            Fronts.Clear();
            int i = 0;
            while (ants.Length != 0)
            {
                List<Ant> front = GetNondominatedIndividuals(ants);
                if (front.Count == 0)
                {
                    Fronts.Add(ants.ToList());
                    break;
                }
                Fronts.Add(front);
                foreach (var path in Fronts[i])
                {
                    ants = ants.Where(val => val != path).ToArray();
                }
                i++;
            }
        }

        public void SystemStart(double distanceTau0, double unbalancingTau0)
        {
            InitializePheromone(distanceTau0, unbalancingTau0);
            
            for (int i = 0; i < _iterations; i++)
            {
                Console.WriteLine("Iteration: " + i);
                InitializeTours();
                
                for (int j = 0; j < _numberOfAnts; j++)
                {
                    TourConstruction(Ants[j]);
                    LocalPheromoneUpdate(Ants[j], distanceTau0, unbalancingTau0);

                    // localSearchResult = LocalSearch(Ants[j].Id, Ants[j].Tour, Ants[j].TourLength);
                    //[j].Tour = localSearchResult.Item1;
                    //[j].TourLength = localSearchResult.Item2;




                    if (Math.Abs(GlobalSecondBestDistance.TourLength) < 0.0000000001 || GlobalSecondBestDistance.TourLength > Ants[j].TourLength)
                    {
                        if (Math.Abs(GlobalBestDistance.TourLength) < 0.0000000001 || GlobalBestDistance.TourLength > Ants[j].TourLength)
                        {
                            GlobalSecondBestDistance = new Ant(GlobalBestDistance);
                            GlobalBestDistance = new Ant(Ants[j]);
                        }
                        else
                        {
                            GlobalSecondBestDistance = new Ant(Ants[j]);
                        }
                        
                    }
                    if (Math.Abs(GlobalSecondBestUnbalancing.UnbalancingDegree) < 0.0000000001 || GlobalSecondBestUnbalancing.UnbalancingDegree > Ants[j].UnbalancingDegree)
                    {
                        if (Math.Abs(GlobalBestUnbalancing.UnbalancingDegree) < 0.0000000001 || GlobalBestUnbalancing.UnbalancingDegree > Ants[j].UnbalancingDegree)
                        {
                            GlobalSecondBestUnbalancing = new Ant(GlobalBestUnbalancing);
                            GlobalBestUnbalancing = new Ant(Ants[j]);
                        }
                        else
                        {
                            GlobalSecondBestUnbalancing = new Ant(Ants[j]);
                        }

                    }
                }
                GetFronts(Ants);
                CollectionOfFirstFronts.AddRange(Fronts[0]);

                GlobalPheromoneUpdate(GlobalBestDistance, GlobalBestUnbalancing, GlobalSecondBestDistance, GlobalSecondBestUnbalancing);
            }
        }

        public void Start()
        {
            var normalizedValues = NormalizeValues();
            double distanceTau0 = 1.0 / normalizedValues.Item1;
            double unbalancingTau0 = 1.0 / normalizedValues.Item2;

            for (int i = 0; i < _numberOfAnts; i++)
            {
                Ants[i] = new Ant(i);
            }
            SystemStart(distanceTau0, unbalancingTau0);
            
        }


    }
}
