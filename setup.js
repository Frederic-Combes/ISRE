function Simulation() {
    this.dimension       = undefined;
    this.time            = {length:2,steps:2000,points:100};
    this.energy          = {count:1000, bins: []};
    this.spinEcho        = {enabled:false, time:[0.25,0.75]};
    this.frequency       = {larmor:0,exchange:0};
    this.inhomogeneity   = {larmor:0,density:0};
    this.relaxation      = 0;
    this.density         = {value:[1], lifetime:Infinity, timeDependence: function (time) {return Math.exp(-time/this.lifetime);}};

    this.name            = undefined;
};

Simulation.prototype.makeList = function (min,max,count)
{
    var array = [];
    for(var val = min; val <= max; val += (max-min)/(count-1))
    {
        array.push(val);
    }
    return array;
}

function c(a,b) { return b.prototype.isPrototypeOf(a);}

Simulation.prototype.check = function()
{
    if( this.dimension !== 1 && this.dimension !== 2 && this.dimension !== 3 )
    {
        print("simulation.dimension must be 1,2 or 3.");
        return false;
    }

    if( typeof this.time.length != "number" || this.time.length <= 0 )
    {
        print("simulation.time.lenght must be a positive number.");
        return false;
    }
    if( typeof this.time.steps != "number" || this.time.steps < 1 )
    {
        print("simulation.time.steps must be a positive number.");
        return false;
    }
    if( typeof this.time.points != "number" || this.time.points < 1 )
    {
        print("simulation.time.points must be a positive number.");
        return false;
    }
    if( this.time.steps % this.time.points != 0)
    {
        print("the number of points to save (simulation.time.points) should divide the number of time steps (simulation.time.steps)");
        this.time.steps = this.time.points * Math.ceil(this.time.steps / this.time.points);
        print("changed the number of time steps (simulation.time.steps) to " + this.time.steps);
    }

    if( typeof this.energy.count != "number" || this.energy.count < 1 )
    {
        print("simulation.energy.count must be a positive number.");
        return false;
    }

    if( c(this.energy.bins, Array) )
    {
        for(var i = 0; i < this.energy.bins.length; i++)
        {
            if( !c(this.energy.bins[i], Array) || this.energy.bins[i].length < 2)
            {
                print("simulation.energy.bins[" + i + "] has to be an array of length 2");
                return false;
            }
        }
    }
    else
    {
        print("simulation.energy.bins has to be an array");
        return false;
    }

    if( typeof this.spinEcho.enabled != "boolean" )
    {
        print("simulation.spinEcho must be either true or false.");
        return false;
    }
    if( this.spinEcho.enabled && typeof this.spinEcho.time != "number" && !c(this.spinEcho.time, Array) )
    {
        print("simulation.spinEcho.time has to be a number or an array when spin echo is enabled.");
        return false;
    }

    if( typeof this.frequency.larmor != "number" && !c(this.frequency.larmor, Array) )
    {
        print("simulation.frequency.larmor is not a number or an array.");
        return false;
    }
    if( typeof this.frequency.exchange != "number" && !c(this.frequency.exchange, Array) )
    {
        print("simulation.frequency.exchange has to be a number or an array.");
        return false;
    }
    if( typeof this.inhomogeneity.larmor != "number" && !c(this.inhomogeneity.larmor, Array) )
    {
        print("simulation.inhomogeneity.larmor has to be a number or an array.");
        return false;
    }
    if( typeof this.inhomogeneity.density != "number" && !c(this.inhomogeneity.density, Array) )
    {
        print("simulation.inhomogeneity.density has to be a number or an array.");
        return false;
    }
    if( typeof this.relaxation != "number" && !c(this.relaxation, Array) )
    {
        print("simulation.relaxation has to be a number or an array.");
        return false;
    }
    if( typeof this.density.value != "number" && !c(this.density.values, Array))
    {
        print("simulation.density.values has to be a number or an array");
        return false;
    }
    if( typeof this.density.timeDependence != "function" )
    {
        print("simulation.density.timeDependence has to be a function of time");
        return false;
    }
    if( typeof this.name != "string" )
    {
        print("simulation.name has to be a file name");
        return false;
    }

    // Turning all single values to arrays
    if( typeof this.frequency.larmor == "number" )
        this.frequency.larmor = [this.frequency.larmor];

    if( typeof this.frequency.exchange == "number" )
        this.frequency.exchange = [this.frequency.exchange];

    if( typeof this.inhomogeneity.larmor == "number" )
        this.inhomogeneity.larmor = [this.inhomogeneity.larmor];

    if( typeof this.inhomogeneity.density == "number" )
        this.inhomogeneity.density = [this.inhomogeneity.density];

    if( typeof this.relaxation == "number" )
        this.relaxation = [this.relaxation];

    if( typeof this.density.value == "number")
        this.density.value = [this.density.value];

    return true;
}

simulation = new Simulation();

simulation.dimension    = 1;

simulation.time.length  = 2;
simulation.time.steps   = 2000;
simulation.time.points  = 100;

simulation.energy.count = 1000;
simulation.energy.bins  = [];

simulation.frequency.larmor         = 0;
simulation.frequency.exchange       = 0;
simulation.inhomogeneity.larmor     = 0;
simulation.inhomogeneity.density    = 0;
simulation.exchange                 = 0;
simulation.relaxation               = 0;

simulation.density.timeDependence   = function(time) {return 1;}




