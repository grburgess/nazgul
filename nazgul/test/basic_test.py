from nazgul import Nazgul



def test_model_building():

    n = Nazgul(assume_effective_area=True)

    n.clean_model()
    
    n = Nazgul(assume_effective_area=False)

    n.clean_model()
    
    n = Nazgul(assume_effective_area=False, earth_occultation=True)

    n.clean_model()
    
    n = Nazgul(assume_effective_area=False, earth_occultation=True, angular_dependence=True)

    n.clean_model()
    
    n = Nazgul(assume_effective_area=False, earth_occultation=False, angular_dependence=True)

    n.clean_model()
    
