function model = load_12ECG_model()

        filename='nine_models.mat';
        A=load(filename);
        model=A.models;

end


