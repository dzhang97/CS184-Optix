from scipy import misc
import numpy as np
from math import floor, ceil

class LightField:
    
    def __init__(self, path_prefix, n, width=480, height=360, bilerp=True):
        self.images = [misc.imread(path_prefix + str(i) + ".png") for i in range(n**2)]
        self.width = width
        self.height = height
        self.n = n
        self.pos = {}
        middle = np.array([int(self.n/2), int(self.n/2)])
        for i in range(n**2):
            self.pos[i] = np.array([i % self.n, i / self.n]) - middle
        self.bilerp = bilerp

    def get_square(self):
        square = []
        for h in range(self.n * self.height):
            row = []
            for w in range(self.n * self.width):
                i = (w % self.n) + self.n * (h % self.n)
                row.append(self.get_value(i, w/self.n, h/self.n))
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

    def get_value(self, i, x, y):
        if (x >= self.width or y >= self.height):
            return np.array([0,0,0])
        return self.images[i][int(y)][int(x)]

    def sample(self, x, y, f):
        color = np.array([0,0,0])
        for i in range(self.n ** 2):
            pos = self.pos[i]
            x_i = x + pos[0]*f
            y_i = y + pos[1]*f
            color += self.get_value(i, x_i, y_i)
        return color / (self.n ** 2)

    def get_refocused(self, f):
        image = []
        for h in range(self.height):
            row = []
            for w in range(self.width):
                row.append(self.sample(w, h, f))
            image.append(row)
        return image

    def shift_image(self, i, x, y):
        image = self.images[i]
        temp = np.zeros_like(image, dtype=np.float64)
        if x < 0:
            temp[:,:x] = image[:,-x:]
        elif x > 0:
            temp[:,x:] = image[:,:-x]
        else:
            temp = image
        ret = np.zeros_like(temp, dtype=np.float64)
        if y < 0:
            ret[:y] = temp[-y:]
        elif y > 0:
            ret[y:] = temp[:-y]
        else:
            ret = temp
        return ret

    def bilerp_shift(self, pos, shiftx, shifty):
        shift_x_int = int(shiftx)
        shift_y_int = int(shifty)
        rx = shiftx - shift_x_int
        ry = shifty - shift_y_int
        if shiftx < 0:
            rx = 1 + rx
        if shifty < 0:
            ry = 1 + ry
        if shiftx >= 0:
            c0 = (1-rx) * self.shift_image(pos, shift_x_int, shift_y_int) \
                + rx * self.shift_image(pos, shift_x_int + 1, shift_y_int)
            c1 = (1-rx) * self.shift_image(pos, shift_x_int, shift_y_int + 1) \
                + rx * self.shift_image(pos, shift_x_int + 1, shift_y_int + 1)
        else:
            c0 = (1-rx) * self.shift_image(pos, shift_x_int - 1, shift_y_int - 1) \
                + rx * self.shift_image(pos, shift_x_int, shift_y_int - 1)
            c1 = (1-rx) * self.shift_image(pos, shift_x_int - 1, shift_y_int) \
                + rx * self.shift_image(pos, shift_x_int, shift_y_int)
        return (1-ry)*c0 + ry*c1

    def refocus_no_skip(self, f, size):
        image = np.zeros(self.images[0].shape)
        shifty = -(f*(self.n-1)/2)
        for j in range(self.n):
            if j > floor((self.n - size) / 2) - 1 and j < self.n - ceil((self.n - size) / 2):
                shiftx = -(f*(self.n-1)/2)
                for i in range(self.n):
                    pos = j*self.n + i
                    if i > floor((self.n - size) / 2) - 1 and i < self.n - ceil((self.n - size) / 2):
                        if self.bilerp:
                            image += self.bilerp_shift(pos, shiftx, shifty)
                        else:
                            image += self.shift_image(pos, int(shiftx), int(shifty))
                    shiftx += f
            shifty += f
        return np.round(image / size**2).astype(np.uint8)

    def refocus_skip(self, f, size, skip):
        old_images, old_n = self.images, self.n
        temp_images = []
        for j in range(int((self.n - skip) // 2), int(self.n - ((self.n - skip) // 2)), int((skip - 1) // 2)):
            for i in range(int((self.n - skip) // 2), int(self.n - ((self.n - skip) // 2)), int((skip - 1) // 2)):
                temp_images.append(self.images[j*self.n + i])
        self.images = temp_images
        self.n = 3
        image = self.refocus_no_skip(f*(skip//2), 3)
        self.images = old_images
        self.n = old_n
        return image

    def refocus(self, f, size, skip):
        if skip == self.n:
            return self.refocus_no_skip(f, size)
        else:
            return self.refocus_skip(f, size, skip)

    def raw_data(self):
        shape = self.images[0].shape
        image = np.zeros((shape[0] * self.n, shape[1] * self.n, shape[2]), dtype=np.uint8)
        for j in range(self.n):
            for i in range(self.n):
                image[i::self.n,j::self.n] = self.images[j*self.n+i]
        return image

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('path')
    parser.add_argument('num_images')
    parser.add_argument('width')
    parser.add_argument('height')
    parser.add_argument('--raw', action='store_true')
    args = parser.parse_args()
    lf = LightField(args.path, int(args.num_images), int(args.width), int(args.height))

    from matplotlib import pyplot as plt
    from matplotlib.widgets import Slider

    if args.raw:
        image = lf.raw_data()
        plt.imshow(image)
        # for val in np.arange(5.0, 11.0, 0.2):
        #     misc.imsave("images/gif/" + str(np.round(val,2)) +".png", lf.refocus(val, 9, 9))

    else:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        fig.subplots_adjust(bottom=.2)
        size = Slider(fig.add_axes([0.15, .09, .75, .03]), 'Num Images', 1, lf.n, valinit=int(lf.n/2), valstep=1)
        shift = Slider(fig.add_axes([0.15, .05, .75, .03]), 'Pixel Shift', 4.0, 12.0, valinit=10.0, valstep=.1)
        skip = Slider(fig.add_axes([0.15, .02, .75, .03]), 'Image Skip', 3, lf.n, valinit=lf.n, valstep=2)
        ax.imshow(lf.refocus(size.val, shift.val, skip.val))
        def update_size(val):
            ax.imshow(lf.refocus(shift.val, val, skip.val))
        def update_focus(val):
            ax.imshow(lf.refocus(val, size.val, skip.val))
        def update_skip(val):
            ax.imshow(lf.refocus(shift.val, size.val, val))
        size.on_changed(update_size)
        shift.on_changed(update_focus)
        skip.on_changed(update_skip)
    plt.show()