import parametry
import seeds

if __name__ == '__main__':
    obj = parametry.mesh("pillar_demo", ext="obj", cls=1, dnt=0)
    obj.emboy(seeds.gen_pillar())
    obj.export()