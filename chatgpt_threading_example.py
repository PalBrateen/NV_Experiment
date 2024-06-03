import threading
import time

# Dummy function to simulate camera behavior
def camera_acquire_frame(trigger_event, num_frames):
    frames_acquired = 0
    while frames_acquired < num_frames:
        trigger_event.wait()  # Wait for the trigger event to be set
        print(f"Camera: Acquiring frame {frames_acquired + 1}")
        frames_acquired += 1
        trigger_event.clear()  # Reset the event to wait for the next pulse

# Main function to handle pulse generation and camera acquisition
def main():
    num_frames = 10
    pulse_interval = 0.5  # Adjust as needed

    # Event object to signal between threads
    trigger_event = threading.Event()

    # Create and start the camera thread
    camera_thread = threading.Thread(target=camera_acquire_frame, args=(trigger_event, num_frames))
    camera_thread.start()

    # Main thread handling pulse generation
    for pulse in range(num_frames):
        time.sleep(pulse_interval)  # Simulate time between pulses
        print(f"Pulse Generator (Main Thread): Sending pulse {pulse + 1}")
        trigger_event.set()  # Set the event to trigger the camera
        trigger_event.clear()

    # Wait for the camera thread to complete
    camera_thread.join()

    print("Acquisition complete")

if __name__ == "__main__":
    main()