using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Collections;

namespace Calculations
{
    //Structs

    //Classes
    public class CSpireMeas
    {
        //Variables
        public bool MIsValid;
        public double MMeasVal;

        //Constructor
        public CSpireMeas()
        {
            MIsValid = false;
            MMeasVal = 0.0;
        }

        //Functions
        //Do this later after Serialization function is done
        public void Serialize()
        {
    
        }

        //missing void as parameter, dont know if necessary yet
        public void InitMeas()
        {
            MIsValid = false;
            MMeasVal = 0.0;
        }

    }
}


