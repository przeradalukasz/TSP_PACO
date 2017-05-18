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
        public Town[] Tour { get; set; }
        public double TourLength { get; set; }

        public Ant(int id)
        {
            Id = id;
        }
    }
}
