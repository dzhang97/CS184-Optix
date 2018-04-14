from scipy import misc

class LightField:
    
    def __init__(self, path_prefix, n, width=480, height=360):
        self.images = [misc.imread(path_prefix + str(i) + ".png") for i in range(n**2)]
        self.width = width
        self.height = height
        self.n = n

    def get_square(self):
        square = []
        for h in range(self.n * self.height):
            row = []
            for w in range(self.n * self.width):
                i = (w % self.n) + self.n * (h % self.n)
                row.append(self.images[i][int(h/self.n)][int(w/self.n)])
            square.append(row)
        return square
        
    def get_side_by_side(self):
        side_by_side = []
        for h in range(self.height):
            row = []
            for w in range(self.width):
                for i in range(self.n**2):
                    row.append(self.images[i][h][w])
            side_by_side.append(row)
        return side_by_side


if __name__ == '__main__':
    from matplotlib import pyplot as plt
    lf = LightField("images/example_image_", 3, 945, 705)
    plt.imshow(lf.get_square())
    plt.show()