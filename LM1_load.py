class Lm1:
    NATIONAL_COEFF = {
        'cl1' : (1.33, 2.40, 1.20, 1.20),
        'cl2' : (1.00, 1.00, 1.00, 1.00)
    }
    def __init__(self, name, load_class, road_width, span_lengths, active_spans=1):
        self.name = name
        self.load_class = load_class
        self.road_width = road_width
        self.span_lengths = span_lengths
        self.active_spans = active_spans
        if self.load_class == 1:
            self.coeffs = Lm1.NATIONAL_COEFF["cl1"]
        else:
            self.coeffs = Lm1.NATIONAL_COEFF["cl2"]
    
    def __str__(self):
        return f'{self.name}'
    
    def data_query(self):
        pass
    
    def compute_total_reactions(self):
        n, w1, l = self.compute_LM1_layout()
        active_spans = list(self.active_spans)
        active_spans = tuple([x-1 for x in active_spans])
        active_length = 0.0
        for i in active_spans:
            active_length += self.span_lengths[i]
        if n == 1:
            lm1_q_reaction = self.coeffs[0] * 9 * active_length * w1 \
                + self.coeffs[2] * 2.5 * active_length * l
        elif n == 2:
            lm1_q_reaction = self.coeffs[0] * 9 * active_length * w1 \
                + self.coeffs[1] * 2.5 * active_length * w1 \
                    + self.coeffs[2] * 2.5 * active_length * l
        else:
            lm1_q_reaction = self.coeffs[0] * 9 * active_length * w1 \
                + self.coeffs[1] * 2.5 * active_length * w1 \
                    + self.coeffs[2] * 2.5 * active_length * (w1+ l)
        ts = [600, 400, 200]
        ts_reaction = sum(ts[0:n])
        total_reaction = lm1_q_reaction + ts_reaction
        return lm1_q_reaction, ts_reaction, total_reaction

    def compute_LM1_layout(self):
        if self.road_width < 5.40:
            numb_lm1_lanes = 1  # n
            width_lm1_lane = 3.00  # w1
            width_leftover_lane = self.road_width - width_lm1_lane  # l
        elif self.road_width < 6.00:
            numb_lm1_lanes = 2
            width_lm1_lane = self.road_width / 2.00
            width_leftover_lane = 0.00
        else:
            numb_lm1_lanes = int(self.road_width // 3)
            width_lm1_lane = 3.00
            width_leftover_lane = self.road_width - numb_lm1_lanes * width_lm1_lane
        return numb_lm1_lanes, width_lm1_lane, width_leftover_lane



def main():
    # nazwa ob, klasa wg EC, szer jezdni, (po kolei dlugosci przesel), (ktore przesla obciazone)
    wd30_s9_lm1 = Lm1('wd30_s9_lm1', 1, 5.7, (28.00, 28.00), (1, 2))
    print(wd30_s9_lm1)
    print(f'liczba pasów, szerokość pasów, pozostałość LM1: {wd30_s9_lm1.compute_LM1_layout()}')
    print(f'reakcja od rozłożonego LM1 {wd30_s9_lm1.compute_total_reactions()[0]}')
    print(f'reakcja od TS LM1 {wd30_s9_lm1.compute_total_reactions()[1]}')
    print(f'reakcja od TS LM1 {wd30_s9_lm1.compute_total_reactions()[2]}')
if __name__ == '__main__':
    main()
