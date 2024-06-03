import asyncio
import cv2
import dcamcon


# Function to capture and display a frame
async def capture_display_frame(cap):
    ret, frame = cap.read()
    if ret:
        cv2.imshow('Live Feed', frame)
        cv2.waitKey(1)  # Wait for 1 millisecond for key press
    else:
        print("Error: Frame not captured")


# Function to handle user input
async def handle_user_input(cap):
    while True:
        # Get user input without blocking the loop
        user_input = await asyncio.get_event_loop().run_in_executor(None, cv2.waitKey)
        # Check for specific key presses (e.g., 'q' to quit)
        if user_input == ord('q'):
            break


# Main function to run everything asynchronously
async def main():
    if dcamcon.dcamcon_init():
        # choose camera and get Dcamcon instance
        dcamcon = dcamcon_choose_and_open()
        if dcamcon is not None:
            res = True
            # set basic properties
            if (not signaled_sigint and
                    res):
                res = setup_properties(dcamcon)

            # show live image
            if (not signaled_sigint and
                    res):
                wait_event_time = show_live_captured_images(dcamcon)

            # close dcam
            dcamcon.close()

    # Open camera capture
    cap = cv2.VideoCapture(0)

    # Create asynchronous tasks
    capture_task = asyncio.create_task(capture_display_frame(cap))
    input_task = asyncio.create_task(handle_user_input(cap))

    # Run tasks concurrently
    await asyncio.gather(capture_task, input_task)

    # Release resources
    cap.release()
    cv2.destroyAllWindows()


if __name__ == '__main__':
    asyncio.run(main())
