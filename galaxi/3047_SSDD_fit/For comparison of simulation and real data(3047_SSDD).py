from PIL import Image, ImageFont, ImageDraw


import time

#filename_real_data = "/home/uliana/Loka/2019/Fitting_3047_SSDD/Real_Data_3047_SSDD_images_other_name/000.jpg" #Real_Data 2949-3021
#filename1 = "/home/uliana/Loka/2019/Fitting_3047_SSDD/2_particles_15nm_Cone/000.jpg" #Cone from 22 layers (15nm, 2 particles)
#filename2 = "/home/uliana/Loka/2019/Fitting_3047_SSDD/2_particles_15nm_Prism/000.jpg" #Prism from 22 layers (15 nm, 2 particles)
#filename3 = "/home/uliana/Loka/2019/Fitting_3047_SSDD/2_particles_mix Prism+Cone_1:1/005.jpg" # Prism+Cone 1:1, 20 layers (15 nm, 2 particles)
#filename4 = "/home/uliana/Loka/2019/Fitting_3047_SSDD/2_particles_mix Prism+Cone_1:3/000.jpg" # Prism+Cone 3:1, 20 layers (15 nm, 2 particles)
#filename5 = "/home/uliana/Loka/2019/Fitting_3047_SSDD/3047_SSDD_Cone6/00.png" # Cone6 only Si, without layers
#filename6 = "/home/uliana/Loka/2019/Fitting_3047_SSDD/3047_SSDD_Cone6+Cap/000.png" # Cone6 (core Si - shall Au structure) and Gold Cap
#filename7 = "/home/uliana/Loka/2019/Fitting_3047_SSDD/3047_SSDD_Prism6/00.png" # Prism6 only Si, without layers
#filename8 = "/home/uliana/Loka/2019/Fitting_3047_SSDD/3047_SSDD_Prism6+Cap/000.png" # Prism6 (core Si - shall Au structure) and Gold Cap
#filename9 = "/home/uliana/Loka/2019/Fitting_3047_SSDD/Cap/1cap_without_tilt/00.jpg" # 1 Gold Cap
#filename10 = "/home/uliana/Loka/2019/Fitting_3047_SSDD/Cap/3caps_tilt30/000.jpg" #3 Gold Caps tilt 30
#filename11 = "/home/uliana/Loka/2019/Fitting_3047_SSDD/Cap/3caps_tilt60/00.jpg" #3 Gold Caps tilt 60
#filename12 = "/home/uliana/Loka/2019/Fitting_3047_SSDD/Prism tilt/3Prism_tilt30_int5+06/000.jpg" #3 silicon prisms, tilt 30
#filename13 = "/home/uliana/Loka/2019/Fitting_3047_SSDD/Prism tilt/3Prism_tilt60_int5+06/000.jpg" #3 silicon prisms, tilt 60



size = 800, 800 #размеры для уменьщения изображения
i=0
while i<120:
    #filename1 = ("/home/uliana/Loka/2019/Fitting_3047_SSDD/2_particles_15nm_Cone/{id:03d}.jpg" .format(id=i)) #Cone from 22 layers (15nm, 2 particles)
    #filename2 = ("/home/uliana/Loka/2019/Fitting_3047_SSDD/2_particles_15nm_Prism/{id:03d}.jpg".format(id=i)) #Prism from 22 layers (15 nm, 2 particles)
    #filename3 = ("/home/uliana/Loka/2019/Fitting_3047_SSDD/2_particles_mix Prism+Cone_1:1/{id:03d}.jpg".format(id=i)) # Prism+Cone 1:1, 20 layers (15 nm, 2 particles)
    #filename4 = ("/home/uliana/Loka/2019/Fitting_3047_SSDD/2_particles_mix Prism+Cone_1:3/{id:03d}.jpg".format(id=i)) # Prism+Cone 3:1, 20 layers (15 nm, 2 particles)
    #filename5 = ("/home/uliana/Loka/2019/Fitting_3047_SSDD/3047_SSDD_Cone6/{id:03d}.png".format(id=i)) # Cone6 only Si, without layers
    #filename6 = ("/home/uliana/Loka/2019/Fitting_3047_SSDD/3047_SSDD_Cone6+Cap/{id:03d}.png".format(id=i)) # Cone6 (core Si - shall Au structure) and Gold Cap
    #filename7 = ("/home/uliana/Loka/2019/Fitting_3047_SSDD/3047_SSDD_Prism6/{id:03d}.png".format(id=i)) # Prism6 only Si, without layers
    #filename8 = ("/home/uliana/Loka/2019/Fitting_3047_SSDD/3047_SSDD_Prism6+Cap/{id:03d}.png".format(id=i)) # Prism6 (core Si - shall Au structure) and Gold Cap
    #filename9 = ("/home/uliana/Loka/2019/Fitting_3047_SSDD/Cap/1cap_without_tilt/000.jpg") # 1 Gold Cap
    #filename10 = ("/home/uliana/Loka/2019/Fitting_3047_SSDD/Cap/3caps_tilt30_other_name/{id:03d}.jpg".format(id=i)) #3 Gold Caps tilt 30
    #filename11 = ("/home/uliana/Loka/2019/Fitting_3047_SSDD/Cap/3caps_tilt60_other_name/{id:03d}.jpg".format(id=i)) #3 Gold Caps tilt 60
    filename12 = ("/home/uliana/Loka/2019/Fitting_3047_SSDD/Prism tilt/3Prism_tilt30_int5+06/{id:03d}.jpg".format(id=i)) #3 silicon prisms, tilt 30
    filename13 = ("/home/uliana/Loka/2019/Fitting_3047_SSDD/Prism tilt/3Prism_tilt60_int5+06/{id:03d}.jpg".format(id=i)) #3 silicon prisms, tilt 60
    filename_real_data = ("/home/uliana/Loka/2019/Fitting_3047_SSDD/Real_Data_3047_SSDD_images_other_name/{id:03d}.jpg" .format(id=i))


    img1 = Image.open(filename_real_data)
    img2 = Image.open(filename12)
    img3 = Image.open(filename13)

    area = (30, 40, 1150, 980)
    img1_cut = img1.crop(area)
    img2_cut = img2.crop(area)
    img3_cut = img3.crop(area)

    img1_cut.thumbnail(size)
    img2_cut.thumbnail(size)
    img3_cut.thumbnail(size)

    img = Image.new('RGB', (800*3, 700))
    img.paste(img1_cut, (0, 0))
    img.paste(img2_cut, (800, 0))
    img.paste(img3_cut, (800*2, 0))
    draw = ImageDraw.Draw(img)
    draw.text(xy=(1200, 20), text = 'Rotation (deg) {id:03d}'.format(id=i), fill=(0, 0, 0))

    img.show()
    i = i + 5
    time.sleep(3)
    #img.save("out.jpg")
    #print("If you would like to see next couple of images, press \"Enter\".")
    #pause = input()


