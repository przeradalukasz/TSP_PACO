using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TSP_PACO
{
    public class Ant
    {
        public int Id { get; set; }
        public List<Town> Tour { get; set; }
        public double TourLength { get; set; }
        public double UnbalancingDegree { get; set; }

        public Ant(Ant ant)
        {
            Id = ant.Id;
            Tour = new List<Town>(ant.Tour);
            TourLength = ant.TourLength;
            UnbalancingDegree = ant.UnbalancingDegree;
        }
        public Ant()
        {
            Tour = new List<Town>();
        }
        public Ant(int id)
        {
            Id = id;
            Tour = new List<Town>();
        }
    }
}
