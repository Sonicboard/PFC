function transducer_object = PFC_TransducerArray2Transducer(transducer_object_list)

transducer_pressure = [];

Tri_Points = [];
Tri_ConnectivityList = [];
Tri_Constraints = [];

for i = 1:numel(transducer_object_list)
    transducer_object_i = transducer_object_list(i);
    
    transducer_pressure_i = transducer_object_i.Pressure;
    
    Tri_Points_i = transducer_object_i.Tri_Points;
    Tri_ConnectivityList_i = transducer_object_i.Tri_ConnectivityList + size(Tri_Points, 1);
    Tri_Constraints_i = transducer_object_i.Tri_Constraints + size(Tri_Points, 1);
    
    transducer_pressure = [transducer_pressure; transducer_pressure_i]; %#ok<AGROW>
    
    Tri_Points = [Tri_Points; Tri_Points_i]; %#ok<AGROW>
    Tri_ConnectivityList = [Tri_ConnectivityList; Tri_ConnectivityList_i]; %#ok<AGROW>
    Tri_Constraints = [Tri_Constraints; Tri_Constraints_i]; %#ok<AGROW>
end

transducer_object.FrequencyDesigned = transducer_object_list(1).FrequencyDesigned;
transducer_object.MediumDesigned = transducer_object_list(1).MediumDesigned;

transducer_object.Pressure = transducer_pressure;

transducer_object.Tri_Points = Tri_Points;
transducer_object.Tri_ConnectivityList = Tri_ConnectivityList;
transducer_object.Tri_Constraints = Tri_Constraints;

end