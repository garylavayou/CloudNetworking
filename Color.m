classdef Color
    properties
      R, G, B
    end

    properties(Dependent)
        RGB
    end
   
    enumeration
        Black       (0,0,0);
        Gray        (0.5,0.5,0.5);
        Blue        (0,0,1);
        Green       (0,1,0);
        Red         (1,0,0);
        Purple      (1,0,1);
        Yellow      (1,1,0);
        White       (1,1,1);
        MildGreen   (0,0.69,0.314);
        MiddleGreen (0.13,0.69,0.3);
        MildBlue    (0,0.45,0.74);
        Orange      (0.85,0.33,0.1);
        DarkPurple	(0.49,0.18,0.56);
        LightGreen  (0.47,0.67,0.19);
        LightBlue	(0.3,0.75,0.93);
    end
    
    methods
        function this = Color(r,g,b)
            if nargin == 1
                this.R = r(1);
                this.G = r(2);
                this.B = r(3);
            else
                this.R = r;
                this.G = g;
                this.B = b;
            end
        end
        
        function rgb = get.RGB(this)
            rgb = [this.R, this.G, this.B];
        end
    end
end