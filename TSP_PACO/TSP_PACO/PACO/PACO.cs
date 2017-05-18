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
        private readonly double _alpha;
        private readonly double _beta;
        private readonly double _q0;
        private readonly int _rho;
        private readonly int _ksi;

        public Town[] Towns;
        private readonly int numberOfCities;
        public double[,] AdjacencyMatrix;
        public double[,] PheromoneMatrix;

        public Ant[] Ants;

        public PACO(double alpha, double beta, double q0, int rho, int ksi, Random rnd)
        {
            _alpha = alpha;
            _beta = beta;
            _q0 = q0;
            _rho = rho;
            _ksi = ksi;
            _rnd = rnd;
        }

    }
}
