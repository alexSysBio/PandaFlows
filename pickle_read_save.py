import pickle

def save_data(data, save_path):
    with open(save_path, 'wb') as handle:
        pickle.dump(data, handle)
        
def load_data(save_path):
    with open(save_path, 'rb') as handle:
        return pickle.load(handle)