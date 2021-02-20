clear
close

% A simple script that converts the NN table from .mat to .txt
% João Inácio, Oct. 7, 2020

type_vals = ["2D_SS", "3D_SC", "3D_BCC", "3D_FCC", "3D_HCP", "3D_Hex"];
NN_vals = [4, 6, 8, 12, 12, 8];

for idx = 1:length(type_vals)
    type = type_vals(idx);
    NN = NN_vals(idx);
    
    for L = 2:8
        load("neighbour_table_" + type + "_1NN_L" + L + ".mat")

        writematrix(NN_table - 1, "txt/neighbour_table_" + type + "_" + NN + "NN_L" + L + ".txt")
        Data = fileread("txt/neighbour_table_" + type + "_" + NN + "NN_L" + L + ".txt");
        Data = strrep(Data, ',', ' ');
        FID = fopen("txt/neighbour_table_" + type + "_" + NN + "NN_L" + L + ".txt", 'w');
        fwrite(FID, Data, 'char');
        fclose(FID);
    end
end
