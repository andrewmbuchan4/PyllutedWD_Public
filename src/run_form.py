print('python3?')
import sys
print(sys.version)
import geology_info as gi
geo_model = gi.GeologyModel()
print(geo_model.form_a_planet_iteratively(66, -2))
