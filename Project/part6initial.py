from AST1100ShortcutSystem import AST1100SolarSystem

seed = 4252
sys = AST1100SolarSystem(seed,hasMoons = True)
sys.landOnPlanet(6,'lander.txt')
