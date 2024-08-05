function medium_data = PFC_GetMedium(medium_name, target_frequency)

medium_data = [];
medium_data.name = medium_name;

if(nargin == 1)
    target_frequency = 40000;  % Initial is 1e6
end

switch(medium_name)
    case 'air'
        medium_data.c = 340;
        medium_data.rho = 1.2;
    case 'water'
        medium_data.c = 1497;
        medium_data.rho = 997;
        
    case 'water_simple'
        medium_data.c = 1500;
        medium_data.rho = 1000;
        
    case 'silicone oil'
        medium_data.c = 980;
        medium_data.rho = 960;
        
    case 'water_loss'
        medium_data.c = 1497;
        medium_data.rho = 997;
        
        medium_data.a = 22e-14; % a (Nepers*s^2/m)
        medium_data.k = 2; % k is dimensionless: loss is exp(-alfa*f^k*z)
        medium_data.alpha = medium_data.a.*target_frequency.^medium_data.k;
        
    case 'soft_tissue'
        medium_data.c = 1497;
        medium_data.rho = 997;
        
        medium_data.a = 7500e-14; % a (Nepers*s^2/m)
        medium_data.k = 2; % k is dimensionless: loss is exp(-alfa*f^k*z)
        medium_data.alpha = medium_data.a.*target_frequency.^medium_data.k;
end

end