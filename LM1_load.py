

class Lm1:
    def __init__(self, name, load_class, road_width, span_lengths):
        self.name = name
        self.load_class = load_class
        self.road_width = road_width
        self.span_lengths = span_lengths
    
    def __str__(self):
        return f'{self.name}'
    
    def data_query(self):
        pass
    
    def compute_total_reactions(self):
        self._compute_LM1_layout()
        total_reactions = self.road_width * self.span_lengths[0] * 2.5
        return total_reactions

    def _compute_LM1_layout(self):
        pass



def main():
    wd30_s9_lm1 = Lm1('wd30_s9_lm1', 2, 7, [28.0, 28.0])
    print(wd30_s9_lm1)
    print(f'całkowita reakcja od LM1: {wd30_s9_lm1.compute_total_reactions()}')
if __name__ == '__main__':
    main()