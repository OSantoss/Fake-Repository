using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Calculations
{
    public class CIirFilter
    {
        //Variables

        protected Globals.Filter MFilter;
        public int MNumCoef;
        public int MNumHist;

        //Constructor
        public CIirFilter(Globals.Filter pFilterTemplate)
        {
            MFilter.Length = pFilterTemplate.Length;
            //no idea why length has to be unsigned int
            MNumCoef = 4 * (int)MFilter.Length + 1;
            MNumHist = 2 * (int)MFilter.Length;

            MFilter.Coef = new double[MNumCoef];
            MFilter.History = new double[MNumHist];

            for(var i = 0; i< MNumCoef;i++)
            {
                MFilter.Coef[i] = pFilterTemplate.Coef[i];
            }

            for(var i=0;i<MNumHist;i++)
            {
                MFilter.History[i] = 0.0;
            }
        }

        public double DoFilter(double input)
        {
            uint i;
            double hist1Ptr, hist2Ptr;
            double output, newHist;
            int coefIndex = 0, hist1Index = 0;

            output = input * MFilter.Coef[coefIndex++];

            for(i = 0;i<MFilter.Length;i++)
            {
                hist1Ptr = MFilter.History[hist1Index];
                hist2Ptr = MFilter.History[hist1Index + 1];

                output = output - hist1Ptr * MFilter.Coef[coefIndex++];
                newHist = output - hist2Ptr * MFilter.Coef[coefIndex++];

                output = newHist + hist1Ptr * MFilter.Coef[coefIndex++];
                output = output + hist2Ptr * MFilter.Coef[coefIndex++];

                MFilter.History[hist1Index] = newHist;
                MFilter.History[hist1Index + 1] = hist1Ptr;
                hist1Index += (int)MFilter.Length;
            }
        
            return output;
        }
    }
}
